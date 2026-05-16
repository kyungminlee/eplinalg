/*
 * qgemm — kind16 (REAL(KIND=16) / __float128) GEMM overlay.
 *
 * Unblocked Netlib reference algorithm with OpenMP parallelism across
 * the outer j loop (columns of C). Every __float128 op lowers to a
 * libquadmath call (~hundreds of cycles), so the arithmetic dominates
 * and packing / cache blocking are unprofitable. Coarse-grain
 * parallelism over j scales nearly linearly with cores.
 *
 * Four algorithmic forms based on (TRANSA, TRANSB), matching Netlib
 * DGEMM. For real kind16, TRANSA='C' ≡ 'T'.
 *
 * Fortran ABI: lowercase name + trailing underscore.
 */

#include <stddef.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

#define QGEMM_OMP_MIN 32

typedef __float128 T;

static int trans_code(const char *p, size_t len) {
    (void)len;
    char c = (char)toupper((unsigned char)*p);
    return (c == 'C') ? 'T' : c;
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]
#define C_(i, j)  c[(size_t)(j) * ldc + (i)]

void qgemm_(
    const char *transa, const char *transb,
    const int *m_, const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    size_t transa_len, size_t transb_len)
{
    const int M = *m_, N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    if (M <= 0 || N <= 0) return;

    const T zero = 0.0Q, one = 1.0Q;

    /* Quick return when only beta scaling is needed. */
    if (alpha == zero || K == 0) {
        for (int j = 0; j < N; ++j) {
            T *cj = &C_(0, j);
            if (beta == zero)      for (int i = 0; i < M; ++i) cj[i]  = zero;
            else if (beta != one)  for (int i = 0; i < M; ++i) cj[i] *= beta;
        }
        return;
    }

#ifdef _OPENMP
    const int use_omp = (N >= QGEMM_OMP_MIN && blas_omp_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        T *cj = &C_(0, j);

        /* beta scaling of column j. */
        if (beta == zero)      for (int i = 0; i < M; ++i) cj[i]  = zero;
        else if (beta != one)  for (int i = 0; i < M; ++i) cj[i] *= beta;

        if (ta == 'N') {
            /* Rank-1 (axpy) form: TEMP = alpha · op(B)[k,j], then
             * C[:,j] += TEMP · A[:,k]. */
            if (tb == 'N') {
                for (int k = 0; k < K; ++k) {
                    const T bkj = B_(k, j);
                    if (bkj != zero) {
                        const T t = alpha * bkj;
                        const T *ak = &A_(0, k);
                        for (int i = 0; i < M; ++i) cj[i] += t * ak[i];
                    }
                }
            } else {  /* TB = 'T' */
                for (int k = 0; k < K; ++k) {
                    const T bjk = B_(j, k);
                    if (bjk != zero) {
                        const T t = alpha * bjk;
                        const T *ak = &A_(0, k);
                        for (int i = 0; i < M; ++i) cj[i] += t * ak[i];
                    }
                }
            }
        } else {  /* TA = 'T': inner-product (DDOT) form. */
            if (tb == 'N') {
                for (int i = 0; i < M; ++i) {
                    T s = zero;
                    for (int k = 0; k < K; ++k) s += A_(k, i) * B_(k, j);
                    cj[i] += alpha * s;
                }
            } else {  /* TB = 'T' */
                for (int i = 0; i < M; ++i) {
                    T s = zero;
                    for (int k = 0; k < K; ++k) s += A_(k, i) * B_(j, k);
                    cj[i] += alpha * s;
                }
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
