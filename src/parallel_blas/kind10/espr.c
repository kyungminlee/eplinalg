/*
 * espr — kind10 (long double) symmetric packed rank-1 update.
 *   A := alpha*x*x^T + A
 */

#include <stddef.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

#define ESPR_OMP_MIN 64

typedef long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

void espr_(
    const char *uplo,
    const int *n_,
    const T *alpha_,
    const T *restrict x, const int *incx_,
    T *restrict ap,
    size_t uplo_len)
{
    (void)uplo_len;
    const int N = *n_;
    const int incx = *incx_;
    const T alpha = *alpha_;
    const T zero = 0.0L;
    const char UPLO = up(uplo);

    if (N == 0 || alpha == zero) return;

    if (incx == 1) {
        if (UPLO == 'U') {
#ifdef _OPENMP
            const int use_omp = (N >= ESPR_OMP_MIN && blas_omp_max_threads() > 1);
            #pragma omp parallel for if(use_omp) schedule(static)
#endif
            for (int j = 0; j < N; ++j) {
                if (x[j] != zero) {
                    const T tmp = alpha * x[j];
                    /* Pointer-walk: gcc emits one `add` + one `cmp`
                     * per iter (9 insns) instead of counter-and-cmp
                     * (10 insns). Matches reference DSPR codegen. */
                    T *restrict apk  = &ap[(size_t)j * (j + 1) / 2];
                    T *restrict aend = apk + j + 1;
                    const T *restrict xp = x;
                    for (; apk < aend; ++apk, ++xp) *apk += *xp * tmp;
                }
            }
        } else {
#ifdef _OPENMP
            const int use_omp = (N >= ESPR_OMP_MIN && blas_omp_max_threads() > 1);
            #pragma omp parallel for if(use_omp) schedule(static)
#endif
            for (int j = 0; j < N; ++j) {
                if (x[j] != zero) {
                    const T tmp = alpha * x[j];
                    T *restrict apk  = &ap[(size_t)j * N - (size_t)j * (j - 1) / 2];
                    T *restrict aend = apk + (N - j);
                    const T *restrict xp = &x[j];
                    for (; apk < aend; ++apk, ++xp) *apk += *xp * tmp;
                }
            }
        }
    } else {
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        int kk = 0;
        if (UPLO == 'U') {
            int jx = kx;
            for (int j = 0; j < N; ++j) {
                if (x[jx] != zero) {
                    const T tmp = alpha * x[jx];
                    int ix = kx;
                    for (int k = kk; k < kk + j + 1; ++k) {
                        ap[k] += x[ix] * tmp;
                        ix += incx;
                    }
                }
                jx += incx;
                kk += j + 1;
            }
        } else {
            int jx = kx;
            for (int j = 0; j < N; ++j) {
                if (x[jx] != zero) {
                    const T tmp = alpha * x[jx];
                    int ix = jx;
                    for (int k = kk; k < kk + N - j; ++k) {
                        ap[k] += x[ix] * tmp;
                        ix += incx;
                    }
                }
                jx += incx;
                kk += N - j;
            }
        }
    }
}
