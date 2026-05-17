/* erot — kind10 real Givens rotation. */
#include <stddef.h>
typedef long double T;

void erot_(const int *n_, T *x, const int *incx_, T *y, const int *incy_,
           const T *c_, const T *s_)
{
    const int n = *n_, incx = *incx_, incy = *incy_;
    const T c = *c_, s = *s_;
    if (n <= 0) return;
    if (incx == 1 && incy == 1) {
        /* Byte-offset shared-index walk: gcc emits one `add` + one
         * `cmp` per iter (16 insns) instead of two pointer increments
         * (18 insns). Matches gfortran reference DROT codegen. See
         * doc/parallel-blas-optimization-findings: Addendum 7. */
        char *restrict xb = (char *)x;
        char *restrict yb = (char *)y;
        const size_t end = (size_t)n * sizeof(T);
        for (size_t k = 0; k < end; k += sizeof(T)) {
            T *xp = (T *)(xb + k);
            T *yp = (T *)(yb + k);
            T tx = c * (*xp) + s * (*yp);
            *yp  = c * (*yp) - s * (*xp);
            *xp  = tx;
        }
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        int iy = (incy < 0) ? (-n + 1) * incy : 0;
        for (int i = 0; i < n; ++i) {
            T tx = c * x[ix] + s * y[iy];
            y[iy] = c * y[iy] - s * x[ix];
            x[ix] = tx;
            ix += incx; iy += incy;
        }
    }
}
