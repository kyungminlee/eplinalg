/* qdot — kind16 real: returns Σ X·Y. */
#include <quadmath.h>
typedef __float128 T;

T qdot_(const int *n_, const T *x, const int *incx_,
        const T *y, const int *incy_)
{
    const int n = *n_, incx = *incx_, incy = *incy_;
    if (n <= 0) return 0.0Q;
    T s0 = 0.0Q, s1 = 0.0Q;
    if (incx == 1 && incy == 1) {
        int i = 0;
        for (; i + 1 < n; i += 2) { s0 += x[i] * y[i]; s1 += x[i + 1] * y[i + 1]; }
        if (i < n) s0 += x[i] * y[i];
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        int iy = (incy < 0) ? (-n + 1) * incy : 0;
        for (int i = 0; i < n; ++i) { s0 += x[ix] * y[iy]; ix += incx; iy += incy; }
    }
    return s0 + s1;
}
