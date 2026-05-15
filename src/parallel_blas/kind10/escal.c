/* escal — kind10 real: X := α · X. x87 fp80, no SIMD path. */
typedef long double T;

void escal_(const int *n_, const T *alpha_, T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    const T alpha = *alpha_;
    if (n <= 0 || alpha == 1.0L) return;
    if (incx == 1) {
        for (int i = 0; i < n; ++i) x[i] *= alpha;
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        for (int i = 0; i < n; ++i) { x[ix] *= alpha; ix += incx; }
    }
}
