/*
 * esymm — kind10 (REAL(KIND=10) / `long double`) symmetric matrix multiply.
 *
 * Blocked "read A_IK once, use twice" pattern (Netlib SYMM logic
 * lifted to block granularity): per off-diagonal block A_IK we fire
 * two consecutive egemm calls that both read the same A_IK from
 * still-warm L2 cache:
 *
 *   USE 1:  C(K, J) += alpha · A_IK^T · B(I, J)   (symmetric reflection)
 *   USE 2:  C(I, J) += alpha · A_IK   · B(K, J)   (direct row I)
 *
 * x87 long double has no SIMD, so peak per-core throughput is bound by
 * the x87 FPU. We size the macro-block so A_IK + an accumulator fit in
 * the per-core L2 between the two gemm calls (16 B/element → 128×128
 * block ≈ 256 KiB), and choose nb so threads each get a J panel.
 */

#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define ESYMM_OMP_MIN 32

typedef long double T;

static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}
/* nb ≈ N/nt for thread balance, capped at 128 so A_IK fits in L2
 * per core (16-byte long double × 128² = 256 KiB), floored at 64
 * to amortize egemm per-call setup. Override via ESYMM_NB. */
static int symm_nb_pick(int M_or_N, int nt) {
    static int g_nb_override = -1;
    if (g_nb_override < 0) g_nb_override = env_int("ESYMM_NB", 0);
    if (g_nb_override > 0) return g_nb_override;
    int nt_eff = (nt > 0) ? nt : 1;
    int nb = (M_or_N + nt_eff - 1) / nt_eff;
    if (nb > 128) nb = 128;
    if (nb < 64)  nb = 64;
    return nb;
}

extern void egemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const T *alpha,
    const T *a, const int *lda,
    const T *b, const int *ldb,
    const T *beta,
    T *c, const int *ldc,
    size_t transa_len, size_t transb_len);

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]
#define C_(i, j)  c[(size_t)(j) * ldc + (i)]

/* Scalar diagonal-block symm for SIDE='L'. */
static void symm_diag_add_L(int ic, int ib, int jc, int jb, T alpha,
                            const T *a, int lda, const T *b, int ldb,
                            T *c, int ldc, char UPLO)
{
    for (int j = jc; j < jc + jb; ++j) {
        T *cj = c + (size_t)j * ldc;
        const T *bj = b + (size_t)j * ldb;
        if (UPLO == 'L') {
            for (int i = ic; i < ic + ib; ++i) {
                const T temp1 = alpha * bj[i];
                T temp2 = 0.0L;
                for (int k = ic; k < i; ++k) {
                    cj[k]  += temp1 * A_(i, k);
                    temp2  += bj[k] * A_(i, k);
                }
                cj[i] += temp1 * A_(i, i) + alpha * temp2;
            }
        } else {
            for (int i = ic + ib - 1; i >= ic; --i) {
                const T temp1 = alpha * bj[i];
                T temp2 = 0.0L;
                for (int k = i + 1; k < ic + ib; ++k) {
                    cj[k]  += temp1 * A_(i, k);
                    temp2  += bj[k] * A_(i, k);
                }
                cj[i] += temp1 * A_(i, i) + alpha * temp2;
            }
        }
    }
}

/* Scalar diagonal-block symm for SIDE='R'. */
static void symm_diag_add_R(int jc, int jb, int ic, int ib, T alpha,
                            const T *a, int lda, const T *b, int ldb,
                            T *c, int ldc, char UPLO)
{
    for (int j = jc; j < jc + jb; ++j) {
        T *cj = c + (size_t)j * ldc;
        {
            const T t = alpha * A_(j, j);
            for (int i = ic; i < ic + ib; ++i) cj[i] += t * B_(i, j);
        }
        if (UPLO == 'L') {
            for (int k = jc; k < j; ++k) {
                const T t = alpha * A_(j, k);
                if (t != 0.0L) for (int i = ic; i < ic + ib; ++i) cj[i] += t * B_(i, k);
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = alpha * A_(k, j);
                if (t != 0.0L) for (int i = ic; i < ic + ib; ++i) cj[i] += t * B_(i, k);
            }
        } else {
            for (int k = jc; k < j; ++k) {
                const T t = alpha * A_(k, j);
                if (t != 0.0L) for (int i = ic; i < ic + ib; ++i) cj[i] += t * B_(i, k);
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = alpha * A_(j, k);
                if (t != 0.0L) for (int i = ic; i < ic + ib; ++i) cj[i] += t * B_(i, k);
            }
        }
    }
}

void esymm_(
    const char *side, const char *uplo,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    size_t side_len, size_t uplo_len)
{
    (void)side_len; (void)uplo_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const char SIDE = up(side);
    const char UPLO = up(uplo);

    if (M == 0 || N == 0) return;

    const T zero = 0.0L, one = 1.0L;
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};

    if (alpha == zero) {
        if (beta == one) return;
#ifdef _OPENMP
        const int axis = (SIDE == 'L') ? N : M;
        const int use_omp = (axis >= ESYMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            T *cj = c + (size_t)j * ldc;
            if (beta == zero) for (int i = 0; i < M; ++i) cj[i] = zero;
            else              for (int i = 0; i < M; ++i) cj[i] *= beta;
        }
        return;
    }

    int nt = 1;
#ifdef _OPENMP
    nt = omp_get_max_threads();
#endif
    const int nb = symm_nb_pick((SIDE == 'L') ? N : M, nt);

    if (SIDE == 'L') {
#ifdef _OPENMP
        const int use_omp = (N >= ESYMM_OMP_MIN && nt > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;

            for (int j = jc; j < jc + jb; ++j) {
                T *cj = c + (size_t)j * ldc;
                if (beta == zero)      for (int i = 0; i < M; ++i) cj[i]  = zero;
                else if (beta != one)  for (int i = 0; i < M; ++i) cj[i] *= beta;
            }

            for (int ic = 0; ic < M; ic += nb) {
                const int ib = (M - ic < nb) ? (M - ic) : nb;

                if (UPLO == 'L') {
                    for (int kc = 0; kc < ic; kc += nb) {
                        const int kb = (ic - kc < nb) ? (ic - kc) : nb;
                        /* USE 1: C(K, J) += alpha · A_IK^T · B(I, J) */
                        egemm_(TN, NN, &kb, &jb, &ib, &alpha,
                               &A_(ic, kc), &lda, &B_(ic, jc), &ldb, &one,
                               &C_(kc, jc), &ldc, 1, 1);
                        /* USE 2: C(I, J) += alpha · A_IK   · B(K, J) */
                        egemm_(NN, NN, &ib, &jb, &kb, &alpha,
                               &A_(ic, kc), &lda, &B_(kc, jc), &ldb, &one,
                               &C_(ic, jc), &ldc, 1, 1);
                    }
                } else {
                    for (int kc = ic + ib; kc < M; kc += nb) {
                        const int kb = (M - kc < nb) ? (M - kc) : nb;
                        egemm_(TN, NN, &kb, &jb, &ib, &alpha,
                               &A_(ic, kc), &lda, &B_(ic, jc), &ldb, &one,
                               &C_(kc, jc), &ldc, 1, 1);
                        egemm_(NN, NN, &ib, &jb, &kb, &alpha,
                               &A_(ic, kc), &lda, &B_(kc, jc), &ldb, &one,
                               &C_(ic, jc), &ldc, 1, 1);
                    }
                }

                symm_diag_add_L(ic, ib, jc, jb, alpha, a, lda, b, ldb, c, ldc, UPLO);
            }
        }
    } else {  /* SIDE = 'R' */
#ifdef _OPENMP
        const int use_omp = (M >= ESYMM_OMP_MIN && nt > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;

            for (int j = 0; j < N; ++j) {
                T *cj = c + (size_t)j * ldc;
                if (beta == zero)      for (int i = ic; i < ic + ib; ++i) cj[i]  = zero;
                else if (beta != one)  for (int i = ic; i < ic + ib; ++i) cj[i] *= beta;
            }

            for (int jc = 0; jc < N; jc += nb) {
                const int jb = (N - jc < nb) ? (N - jc) : nb;

                if (UPLO == 'L') {
                    for (int kc = jc + jb; kc < N; kc += nb) {
                        const int kb = (N - kc < nb) ? (N - kc) : nb;
                        egemm_(NN, NN, &ib, &jb, &kb, &alpha,
                               &B_(ic, kc), &ldb, &A_(kc, jc), &lda, &one,
                               &C_(ic, jc), &ldc, 1, 1);
                        egemm_(NN, TN, &ib, &kb, &jb, &alpha,
                               &B_(ic, jc), &ldb, &A_(kc, jc), &lda, &one,
                               &C_(ic, kc), &ldc, 1, 1);
                    }
                } else {
                    for (int kc = 0; kc < jc; kc += nb) {
                        const int kb = (jc - kc < nb) ? (jc - kc) : nb;
                        egemm_(NN, NN, &ib, &jb, &kb, &alpha,
                               &B_(ic, kc), &ldb, &A_(kc, jc), &lda, &one,
                               &C_(ic, jc), &ldc, 1, 1);
                        egemm_(NN, TN, &ib, &kb, &jb, &alpha,
                               &B_(ic, jc), &ldb, &A_(kc, jc), &lda, &one,
                               &C_(ic, kc), &ldc, 1, 1);
                    }
                }

                symm_diag_add_R(jc, jb, ic, ib, alpha, a, lda, b, ldb, c, ldc, UPLO);
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
