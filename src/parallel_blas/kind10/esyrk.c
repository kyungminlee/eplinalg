/*
 * esyrk — kind10 (REAL(KIND=10) / `long double`) symmetric rank-k update.
 *
 *   C := alpha · A · Aᵀ + beta · C          (TRANS='N')
 *   C := alpha · Aᵀ · A + beta · C          (TRANS='T'/'C')
 *
 * C is N×N symmetric; only the UPLO triangle is touched.
 *
 * Blocked: each thread owns one nb-wide column block of C and does
 *   (1) beta-scale its UPLO slice,
 *   (2) scalar rank-k on the jb×jb diagonal block,
 *   (3) one egemm call for the trailing (UPLO=L: rows below; UPLO=U:
 *       rows above) — this is where the perf win comes from at kind10,
 *       since egemm is ~2× the unblocked rank-1.
 *
 * Outer omp parallel-for over jc blocks with dynamic scheduling so the
 * imbalanced "trailing height varies with jc" pattern stays smooth
 * across threads.
 */

#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define ESYRK_OMP_MIN 32

typedef long double T;

static int g_esyrk_nb = 64;

__attribute__((constructor))
static void esyrk_init(void) {
    const char *s = getenv("ESYRK_NB");
    if (s && *s) {
        int v = atoi(s);
        if (v > 0) g_esyrk_nb = v;
    }
}
static int syrk_nb(void) { return g_esyrk_nb; }

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
#define C_(i, j)  c[(size_t)(j) * ldc + (i)]

/* Rank-k addition into the UPLO triangle of the diagonal jb×jb block
 * at C[jc..jc+jb, jc..jc+jb]. No beta scaling (caller pre-scales). */
static void syrk_diag_add(int jc, int jb, int K, T alpha,
                          const T *restrict a, int lda,
                          T *restrict c, int ldc,
                          char UPLO, int TR)
{
    if (TR == 'N') {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + (size_t)j * ldc;
            for (int l = 0; l < K; ++l) {
                const T ajl = A_(j, l);
                if (ajl != 0.0L) {
                    const T t = alpha * ajl;
                    for (int i = i_lo; i < i_hi; ++i) cj[i] += t * A_(i, l);
                }
            }
        }
    } else {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + (size_t)j * ldc;
            const T *Aj = a + (size_t)j * lda;
            for (int i = i_lo; i < i_hi; ++i) {
                const T *Ai = a + (size_t)i * lda;
                T s = 0.0L;
                for (int l = 0; l < K; ++l) s += Ai[l] * Aj[l];
                cj[i] += alpha * s;
            }
        }
    }
}

void esyrk_(
    const char *uplo, const char *trans,
    const int *n_, const int *k_,
    const T *alpha_,
    const T *restrict a, const int *lda_,
    const T *beta_,
    T *restrict c, const int *ldc_,
    size_t uplo_len, size_t trans_len)
{
    (void)uplo_len; (void)trans_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);
    char TR = up(trans);
    if (TR == 'C') TR = 'T';

    if (N == 0) return;

    const T zero = 0.0L, one = 1.0L;
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};

    /* alpha or K trivial — only beta-scale. */
    if (alpha == zero || K == 0) {
        if (beta == one) return;
#ifdef _OPENMP
        const int use_omp = (N >= ESYRK_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + (size_t)j * ldc;
            if (beta == zero) for (int i = i_lo; i < i_hi; ++i) cj[i] = zero;
            else              for (int i = i_lo; i < i_hi; ++i) cj[i] *= beta;
        }
        return;
    }

    const int nb = syrk_nb();

#ifdef _OPENMP
    const int use_omp = (N >= ESYRK_OMP_MIN && omp_get_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
    for (int jc = 0; jc < N; jc += nb) {
        const int jb = (N - jc < nb) ? (N - jc) : nb;

        /* (1) Beta-scale this block's UPLO slice of C. */
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + (size_t)j * ldc;
            if (beta == zero)      for (int i = i_lo; i < i_hi; ++i) cj[i]  = zero;
            else if (beta != one)  for (int i = i_lo; i < i_hi; ++i) cj[i] *= beta;
        }

        /* (2) Diagonal block: scalar rank-k add. */
        syrk_diag_add(jc, jb, K, alpha, a, lda, c, ldc, UPLO, TR);

        /* (3) Off-diagonal trailing via egemm with beta=1. */
        if (UPLO == 'L') {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int j0 = jc + jb;
                if (TR == 'N') {
                    egemm_(NN, TN, &trailing, &jb, &K, &alpha,
                           &A_(j0, 0), &lda,
                           &A_(jc, 0), &lda,
                           &one,
                           &C_(j0, jc), &ldc, 1, 1);
                } else {
                    egemm_(TN, NN, &trailing, &jb, &K, &alpha,
                           &A_(0, j0), &lda,
                           &A_(0, jc), &lda,
                           &one,
                           &C_(j0, jc), &ldc, 1, 1);
                }
            }
        } else {  /* UPLO == 'U' */
            if (jc > 0) {
                if (TR == 'N') {
                    egemm_(NN, TN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda,
                           &A_(jc, 0), &lda,
                           &one,
                           &C_(0, jc), &ldc, 1, 1);
                } else {
                    egemm_(TN, NN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda,
                           &A_(0, jc), &lda,
                           &one,
                           &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef C_
