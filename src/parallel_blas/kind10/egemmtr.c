/*
 * egemmtr — kind10 real (long double) triangular GEMM update.
 *   C := alpha · op(A) · op(B) + beta · C   (only UPLO triangle of C)
 *
 * Structure mirrors esyrk / mgemmtr: blocked over jc.
 *   - jb×jb diagonal triangle: scalar rank-k add.
 *   - Off-diagonal rectangle: routed through egemm, which carries
 *     the GotoBLAS-style packing + MR×NR x87 register tiles. This is
 *     where most of the work lives once N ≳ 3·nb; without this hand-off
 *     the unblocked loop would throw away egemm's packing win.
 * Dynamic outer schedule because per-block trailing height varies
 * with jc (same imbalance pattern as esyrk).
 */
#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define EGEMMTR_OMP_MIN 32

typedef long double T;

static int g_egemmtr_nb = 0;

static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}

static int gemmtr_nb(void) {
    if (g_egemmtr_nb == 0) g_egemmtr_nb = env_int("EGEMMTR_NB", 64);
    return g_egemmtr_nb;
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
static int trans_code(const char *p) {
    char c = up(p);
    return (c == 'C') ? 'T' : c;
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]
#define C_(i, j)  c[(size_t)(j) * ldc + (i)]

/* Scalar update of the jb×jb diagonal triangle at (jc, jc).
 * Assumes beta-scaling on C[is..ie, j] already done. */
static void diag_add(int jc, int jb, int K, T alpha,
                     const T *restrict a, int lda,
                     const T *restrict b, int ldb,
                     T *restrict c, int ldc,
                     int upper, char ta, char tb)
{
    const T zero = 0.0L;
    for (int j = jc; j < jc + jb; ++j) {
        const int is = upper ? jc        : j;
        const int ie = upper ? (j + 1)   : (jc + jb);
        T *cj = c + (size_t)j * ldc;

        if (ta == 'N') {
            if (tb == 'N') {
                for (int l = 0; l < K; ++l) {
                    const T t = alpha * B_(l, j);
                    if (t != zero) {
                        const T *al = &A_(0, l);
                        for (int i = is; i < ie; ++i) cj[i] += t * al[i];
                    }
                }
            } else {
                for (int l = 0; l < K; ++l) {
                    const T t = alpha * B_(j, l);
                    if (t != zero) {
                        const T *al = &A_(0, l);
                        for (int i = is; i < ie; ++i) cj[i] += t * al[i];
                    }
                }
            }
        } else {
            if (tb == 'N') {
                for (int i = is; i < ie; ++i) {
                    T s = zero;
                    for (int l = 0; l < K; ++l) s += A_(l, i) * B_(l, j);
                    cj[i] += alpha * s;
                }
            } else {
                for (int i = is; i < ie; ++i) {
                    T s = zero;
                    for (int l = 0; l < K; ++l) s += A_(l, i) * B_(j, l);
                    cj[i] += alpha * s;
                }
            }
        }
    }
}

void egemmtr_(const char *uplo, const char *transa, const char *transb,
              const int *n_, const int *k_,
              const T *alpha_,
              const T *a, const int *lda_,
              const T *b, const int *ldb_,
              const T *beta_,
              T *c, const int *ldc_,
              size_t uplo_len, size_t ta_len, size_t tb_len)
{
    (void)uplo_len; (void)ta_len; (void)tb_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int upper = (up(uplo) == 'U');
    const char ta = (char)trans_code(transa);
    const char tb = (char)trans_code(transb);

    if (N <= 0) return;
    const T zero = 0.0L, one = 1.0L;

    if (alpha == zero || K == 0) {
        if (beta == one) return;
#ifdef _OPENMP
        const int use_omp0 = (N >= EGEMMTR_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp0) schedule(static, 1)
#endif
        for (int j = 0; j < N; ++j) {
            const int is = upper ? 0 : j;
            const int ie = upper ? (j + 1) : N;
            T *cj = &C_(0, j);
            if (beta == zero) for (int i = is; i < ie; ++i) cj[i]  = zero;
            else              for (int i = is; i < ie; ++i) cj[i] *= beta;
        }
        return;
    }

    const int nb = gemmtr_nb();
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    const char *ta_s = (ta == 'N') ? NN : TN;
    const char *tb_s = (tb == 'N') ? NN : TN;

#ifdef _OPENMP
    const int use_omp = (N >= EGEMMTR_OMP_MIN && omp_get_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
    for (int jc = 0; jc < N; jc += nb) {
        const int jb = (N - jc < nb) ? (N - jc) : nb;

        /* (1) Beta-scale the triangle slice for cols [jc, jc+jb). */
        for (int j = jc; j < jc + jb; ++j) {
            const int is = upper ? 0 : j;
            const int ie = upper ? (j + 1) : N;
            T *cj = &C_(0, j);
            if (beta == zero)      for (int i = is; i < ie; ++i) cj[i]  = zero;
            else if (beta != one)  for (int i = is; i < ie; ++i) cj[i] *= beta;
        }

        /* (2) Diagonal jb×jb triangle: scalar. */
        diag_add(jc, jb, K, alpha, a, lda, b, ldb, c, ldc, upper, ta, tb);

        /* (3) Off-diagonal rectangle: routed through egemm (packing). */
        if (upper) {
            if (jc > 0) {
                const int m = jc;
                const T *ablk = (ta == 'N') ? &A_(0, 0) : &A_(0, 0);
                const T *bblk = (tb == 'N') ? &B_(0, jc) : &B_(jc, 0);
                egemm_(ta_s, tb_s, &m, &jb, &K, &alpha,
                       ablk, &lda, bblk, &ldb,
                       &one, &C_(0, jc), &ldc, 1, 1);
            }
        } else {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int r0 = jc + jb;
                const T *ablk = (ta == 'N') ? &A_(r0, 0) : &A_(0, r0);
                const T *bblk = (tb == 'N') ? &B_(0, jc) : &B_(jc, 0);
                egemm_(ta_s, tb_s, &trailing, &jb, &K, &alpha,
                       ablk, &lda, bblk, &ldb,
                       &one, &C_(r0, jc), &ldc, 1, 1);
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
