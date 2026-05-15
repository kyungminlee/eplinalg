/* yescal — kind10: X := α·X with α real, X complex (CSSCAL/ZDSCAL analog). */
typedef _Complex long double T;
typedef long double R;

void yescal_(const int *n_, const R *alpha_, T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    const R alpha = *alpha_;
    if (n <= 0 || alpha == 1.0L) return;
    if (incx == 1) {
        for (int i = 0; i < n; ++i)
            x[i] = (__real__ x[i] * alpha) + (__imag__ x[i] * alpha) * 1.0iL;
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        for (int i = 0; i < n; ++i) {
            x[ix] = (__real__ x[ix] * alpha) + (__imag__ x[ix] * alpha) * 1.0iL;
            ix += incx;
        }
    }
}
