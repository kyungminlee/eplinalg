/*
 * xgemm — kind16 complex (__complex128 / COMPLEX(KIND=16)) GEMM overlay.
 *
 * Unblocked Netlib reference (ZGEMM-style) with OpenMP parallelism
 * across columns of C. No packing or cache blocking — every op is a
 * libquadmath call, the arithmetic dominates, and coarse-grain
 * parallelism scales nearly linearly with cores.
 *
 * TRANSA / TRANSB independently in {N, T, C}. Conjugation under 'C'
 * is applied at element access time; the per-element branch is one
 * sign-flip on the imag part, dwarfed by the libquadmath cost of the
 * surrounding multiply, so the conj flag costs effectively nothing.
 */

#include <stddef.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

#define XGEMM_OMP_MIN 32

typedef __complex128 T;

static int trans_code(const char *p, size_t len) {
    (void)len;
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]
#define C_(i, j)  c[(size_t)(j) * ldc + (i)]

void xgemm_(
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

    const T zero = 0.0Q + 0.0Qi;
    const T one  = 1.0Q + 0.0Qi;

    if (alpha == zero || K == 0) {
        for (int j = 0; j < N; ++j) {
            T *cj = &C_(0, j);
            if (beta == zero)      for (int i = 0; i < M; ++i) cj[i]  = zero;
            else if (beta != one)  for (int i = 0; i < M; ++i) cj[i] *= beta;
        }
        return;
    }

    const int conj_a = (ta == 'C');
    const int conj_b = (tb == 'C');
    const int trans_a = (ta != 'N');
    const int trans_b = (tb != 'N');

#ifdef _OPENMP
    const int use_omp = (N >= XGEMM_OMP_MIN && blas_omp_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        T *cj = &C_(0, j);

        if (beta == zero)      for (int i = 0; i < M; ++i) cj[i]  = zero;
        else if (beta != one)  for (int i = 0; i < M; ++i) cj[i] *= beta;

        if (!trans_a) {
            /* op(A) = A: axpy form. */
            for (int k = 0; k < K; ++k) {
                T bkj;
                if (!trans_b)      bkj = B_(k, j);
                else if (!conj_b)  bkj = B_(j, k);
                else               bkj = conjq(B_(j, k));
                if (bkj != zero) {
                    const T t = alpha * bkj;
                    const T *ak = &A_(0, k);
                    for (int i = 0; i < M; ++i) cj[i] += t * ak[i];
                }
            }
        } else {
            /* op(A) ∈ {Aᵀ, Aᴴ}: inner-product form. */
            for (int i = 0; i < M; ++i) {
                T s = zero;
                if (!trans_b) {
                    if (!conj_a) for (int k = 0; k < K; ++k) s += A_(k, i)         * B_(k, j);
                    else         for (int k = 0; k < K; ++k) s += conjq(A_(k, i))  * B_(k, j);
                } else if (!conj_b) {
                    if (!conj_a) for (int k = 0; k < K; ++k) s += A_(k, i)         * B_(j, k);
                    else         for (int k = 0; k < K; ++k) s += conjq(A_(k, i))  * B_(j, k);
                } else {
                    if (!conj_a) for (int k = 0; k < K; ++k) s += A_(k, i)         * conjq(B_(j, k));
                    else         for (int k = 0; k < K; ++k) s += conjq(A_(k, i))  * conjq(B_(j, k));
                }
                cj[i] += alpha * s;
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
