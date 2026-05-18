/*
 * eger — kind10 (REAL(KIND=10)) general rank-1 update.
 *   A := alpha · x · yᵀ + A    where A is M×N
 *
 * Bandwidth-bound (~1 flop per A element). The structural payoff is
 * OMP-over-j (each column independent) + restrict-based alias info.
 */

#include <stddef.h>
#include <ctype.h>
#include "../common/blas_omp.h"

#define EGER_OMP_MIN 64

typedef long double T;

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void eger_(
    const int *m_, const int *n_,
    const T *alpha_,
    const T *restrict x, const int *incx_,
    const T *restrict y, const int *incy_,
    T *restrict a, const int *lda_)
{
    const int M = *m_, N = *n_;
    const int incx = *incx_, incy = *incy_, lda = *lda_;
    const T alpha = *alpha_;
    const T zero = 0.0L;

    if (M == 0 || N == 0 || alpha == zero) return;

    if (incx == 1 && incy == 1) {
        /* blas_omp_max_threads is cached — avoids paying the libgomp
         * function call on every BLAS invocation (was a ~15% hit at
         * small N for x87 long double work). */
        const int use_omp = (N >= EGER_OMP_MIN && blas_omp_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
        for (int j = 0; j < N; ++j) {
            const T yj = y[j];
            if (yj != zero) {
                const T t = alpha * yj;
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) aj[i] += t * x[i];
            }
        }
    } else {
        int jy = (incy < 0) ? -(N - 1) * incy : 0;
        for (int j = 0; j < N; ++j) {
            const T yj = y[jy];
            if (yj != zero) {
                const T t = alpha * yj;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) {
                    aj[i] += t * x[ix];
                    ix += incx;
                }
            }
            jy += incy;
        }
    }
}

#undef A_
