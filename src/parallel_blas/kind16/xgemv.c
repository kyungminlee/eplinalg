/*
 * xgemv — kind16 complex (__complex128) general matrix-vector multiply.
 */

#include <stddef.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

#define XGEMV_OMP_MIN 64

typedef __complex128 T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

static inline T cconj(T z) { return conjq(z); }

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void xgemv_(
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
    const char TR = up(trans);

    if (M == 0 || N == 0) return;

    const T zero = 0.0Q + 0.0Qi, one = 1.0Q + 0.0Qi;
    const int leny = (TR == 'N') ? M : N;

    if (beta != one) {
        int iy = (incy < 0) ? -(leny - 1) * incy : 0;
        for (int i = 0; i < leny; ++i) {
            if (beta == zero) y[iy] = zero;
            else              y[iy] *= beta;
            iy += incy;
        }
    }

    if (alpha == zero) return;

    if (TR == 'N' && incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (M >= XGEMV_OMP_MIN && blas_omp_max_threads() > 1);
        #pragma omp parallel if(use_omp)
        {
            int tid = 0, nt = 1;
            if (use_omp) { tid = omp_get_thread_num(); nt = omp_get_num_threads(); }
            const int i_lo = ((long long)M * tid) / nt;
            const int i_hi = ((long long)M * (tid + 1)) / nt;
            for (int j = 0; j < N; ++j) {
                const T xj = x[j];
                if (xj != zero) {
                    const T t = alpha * xj;
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] += t * aj[i];
                }
            }
        }
#else
        for (int j = 0; j < N; ++j) {
            const T xj = x[j];
            if (xj != zero) {
                const T t = alpha * xj;
                const T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) y[i] += t * aj[i];
            }
        }
#endif
    } else if ((TR == 'T' || TR == 'C') && incx == 1 && incy == 1) {
        const int conj_a = (TR == 'C');
#ifdef _OPENMP
        const int use_omp = (N >= XGEMV_OMP_MIN && blas_omp_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const T *aj = &A_(0, j);
            T s = zero;
            if (conj_a) {
                for (int i = 0; i < M; ++i) s += cconj(aj[i]) * x[i];
            } else {
                for (int i = 0; i < M; ++i) s += aj[i] * x[i];
            }
            y[j] += alpha * s;
        }
    } else {
        if (TR == 'N') {
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
        } else {
            const int conj_a = (TR == 'C');
            int jy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int j = 0; j < N; ++j) {
                T s = zero;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                for (int i = 0; i < M; ++i) {
                    s += (conj_a ? cconj(A_(i, j)) : A_(i, j)) * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * s;
                jy += incy;
            }
        }
    }
}

#undef A_
