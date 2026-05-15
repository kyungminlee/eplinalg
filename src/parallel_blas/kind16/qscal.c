/* qscal — kind16 real: X := α · X. */
#include <quadmath.h>
typedef __float128 T;

void qscal_(const int *n_, const T *alpha_, T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    const T alpha = *alpha_;
    if (n <= 0 || alpha == 1.0Q) return;
    if (incx == 1) for (int i = 0; i < n; ++i) x[i] *= alpha;
    else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        for (int i = 0; i < n; ++i) { x[ix] *= alpha; ix += incx; }
    }
}
