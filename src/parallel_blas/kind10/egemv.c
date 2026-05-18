/*
 * egemv — kind10 (REAL(KIND=10) / `long double`) general matrix-vector.
 *
 *   y := alpha · A · x + beta · y          (TRANS='N')   A is M×N
 *   y := alpha · Aᵀ · x + beta · y         (TRANS='T'/'C') A is M×N, y is N
 *
 * Reference DGEMV algorithm. L2 is bandwidth-bound (~1 flop/element),
 * so structural optimization is limited; the overlay's value is
 * OpenMP parallelism (gfortran's reference DGEMV is serial) and
 * restrict-based alias elimination.
 *
 * Partitioning strategy: each thread owns a contiguous slice of the
 * OUTPUT vector y. For TRANS='N' that's rows of A (and y has length M).
 * For TRANS='T'/'C' that's columns of A (and y has length N). All A
 * accesses inside a thread's slice remain column-major stride-1.
 */

#include <stddef.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

#define EGEMV_OMP_MIN 64

typedef long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void egemv_(
    const char *trans,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *restrict a, const int *lda_,
    const T *restrict x, const int *incx_,
    const T *beta_,
    T *restrict y, const int *incy_,
    size_t trans_len)
{
    (void)trans_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, incx = *incx_, incy = *incy_;
    const T alpha = *alpha_, beta = *beta_;
    char TR = up(trans);
    if (TR == 'C') TR = 'T';

    if (M == 0 || N == 0) return;

    const T zero = 0.0L, one = 1.0L;

    /* Output length: M for TRANS='N', N for TRANS='T'. */
    const int leny = (TR == 'N') ? M : N;

    /* Beta-scale y (or zero) first — independent of A·x update. */
    if (beta != one) {
        if (incy == 1) {
            if (beta == zero) for (int i = 0; i < leny; ++i) y[i] = zero;
            else              for (int i = 0; i < leny; ++i) y[i] *= beta;
        } else {
            size_t iy0 = (incy < 0) ? (size_t)(-(leny - 1)) * (size_t)incy : 0;
            (void)iy0;
            int iy = (incy < 0) ? -(leny - 1) * incy : 0;
            for (int i = 0; i < leny; ++i) {
                if (beta == zero) y[iy] = zero;
                else              y[iy] *= beta;
                iy += incy;
            }
        }
    }

    if (alpha == zero) return;

    if (TR == 'N') {
        /* y[i] += alpha · sum_j A(i,j) · x(j).
         * Stride-1 contiguous case (most common: incx==incy==1).
         * Partition rows of A across threads so each thread updates
         * its own slice of y with no cross-thread writes. */
        if (incx == 1 && incy == 1) {
            /* J-axis unroll by 2: process two columns at a time so each
             * y[i] load+store services both column updates. Matches
             * gfortran's reference DGEMV codegen — halves y memory
             * traffic on the column-AXPY path. */
#ifdef _OPENMP
            const int use_omp = (M >= EGEMV_OMP_MIN && blas_omp_max_threads() > 1);
#endif
            int i_lo = 0, i_hi = M;
            (void)i_lo; (void)i_hi;
#ifdef _OPENMP
            #pragma omp parallel if(use_omp) firstprivate(i_lo, i_hi)
            {
                if (use_omp) {
                    const int tid = omp_get_thread_num();
                    const int nt  = omp_get_num_threads();
                    i_lo = ((long long)M * tid) / nt;
                    i_hi = ((long long)M * (tid + 1)) / nt;
                }
#endif
                const int span = i_hi - i_lo;
                int j = 0;
                for (; j + 1 < N; j += 2) {
                    const T t0 = alpha * x[j];
                    const T t1 = alpha * x[j + 1];
                    /* Shared-index walk: one base pointer + a byte
                     * offset that's added once per iter, used as the
                     * `index` operand for all three loads. Mirrors
                     * gfortran reference DGEMV codegen (one `add` +
                     * one `cmp` per iter; three bases share a single
                     * index register). Parenthesization on the inner
                     * expression matches gfortran's interleaved
                     * fmul/fadd schedule. */
                    /* gcc emits a shared-index pointer-walk only when
                     * the inner loop is structured as one byte offset
                     * `k` added to three different base pointers. With
                     * the natural `for (i ...) yp[i] = ... a0[i]
                     * ... a1[i]` shape, gcc keeps three independent
                     * pointers and emits 3 `add`s per iter (13 insns).
                     * The byte-offset form drops to 11 insns per iter
                     * and matches gfortran reference DGEMV codegen. */
                    char *restrict yp = (char *)(y + i_lo);
                    const char *restrict a0 = (const char *)&A_(i_lo, j);
                    const char *restrict a1 = (const char *)&A_(i_lo, j + 1);
                    const size_t end = (size_t)span * sizeof(T);
                    for (size_t k = 0; k < end; k += sizeof(T)) {
                        T *yk        = (T *)(yp + k);
                        const T *a0k = (const T *)(a0 + k);
                        const T *a1k = (const T *)(a1 + k);
                        *yk = (*yk + t0 * *a0k) + t1 * *a1k;
                    }
                }
                for (; j < N; ++j) {
                    const T t = alpha * x[j];
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] += t * aj[i];
                }
#ifdef _OPENMP
            }
#endif
        } else if (incy == 1) {
            /* incx != 1, incy == 1: y access is unit-stride, so the same
             * J-axis unroll-by-2 shared-index walk applies — just read
             * x[jx] with stride. Matches Fortran reference's (INCY.EQ.1)
             * inner branch. Without this case the strided-x cells
             * (N/x-1, N/x2) collapsed to ~0.64× of migrated because the
             * old "general stride" fallback used `iy += incy` even when
             * incy==1, preventing J-unroll. */
#ifdef _OPENMP
            const int use_omp = (M >= EGEMV_OMP_MIN && blas_omp_max_threads() > 1);
#endif
            int i_lo = 0, i_hi = M;
            (void)i_lo; (void)i_hi;
#ifdef _OPENMP
            #pragma omp parallel if(use_omp) firstprivate(i_lo, i_hi)
            {
                if (use_omp) {
                    const int tid = omp_get_thread_num();
                    const int nt  = omp_get_num_threads();
                    i_lo = ((long long)M * tid) / nt;
                    i_hi = ((long long)M * (tid + 1)) / nt;
                }
#endif
                const int span = i_hi - i_lo;
                int jx = (incx < 0) ? -(N - 1) * incx : 0;
                int j = 0;
                for (; j + 1 < N; j += 2) {
                    const T t0 = alpha * x[jx];
                    const T t1 = alpha * x[jx + incx];
                    char *restrict yp = (char *)(y + i_lo);
                    const char *restrict a0 = (const char *)&A_(i_lo, j);
                    const char *restrict a1 = (const char *)&A_(i_lo, j + 1);
                    const size_t end = (size_t)span * sizeof(T);
                    for (size_t k = 0; k < end; k += sizeof(T)) {
                        T *yk        = (T *)(yp + k);
                        const T *a0k = (const T *)(a0 + k);
                        const T *a1k = (const T *)(a1 + k);
                        *yk = (*yk + t0 * *a0k) + t1 * *a1k;
                    }
                    jx += 2 * incx;
                }
                for (; j < N; ++j) {
                    const T t = alpha * x[jx];
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] += t * aj[i];
                    jx += incx;
                }
#ifdef _OPENMP
            }
#endif
        } else {
            /* incy != 1: strided y writes — naive reference shape. */
            int jx = (incx < 0) ? -(N - 1) * incx : 0;
            for (int j = 0; j < N; ++j) {
                const T xj = x[jx];
                if (xj != zero) {
                    const T t = alpha * xj;
                    int iy = (incy < 0) ? -(M - 1) * incy : 0;
                    for (int i = 0; i < M; ++i) {
                        y[iy] += t * A_(i, j);
                        iy += incy;
                    }
                }
                jx += incx;
            }
        }
    } else {  /* TRANS = 'T' or 'C' (real: same): y[j] += alpha · sum_i A(i,j) · x(i) */
        if (incx == 1 && incy == 1) {
#ifdef _OPENMP
            const int use_omp = (N >= EGEMV_OMP_MIN && blas_omp_max_threads() > 1);
            #pragma omp parallel for if(use_omp) schedule(static)
#endif
            for (int j = 0; j < N; ++j) {
                const T *aj = &A_(0, j);
                /* 2-chain dot product for x87 pipelining. */
                T s0 = zero, s1 = zero;
                int i = 0;
                for (; i + 1 < M; i += 2) {
                    s0 += aj[i]     * x[i];
                    s1 += aj[i + 1] * x[i + 1];
                }
                T s = s0 + s1;
                for (; i < M; ++i) s += aj[i] * x[i];
                y[j] += alpha * s;
            }
        } else {
            int jy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int j = 0; j < N; ++j) {
                T s = zero;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                for (int i = 0; i < M; ++i) {
                    s += A_(i, j) * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * s;
                jy += incy;
            }
        }
    }
}

#undef A_
