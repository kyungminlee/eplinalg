/* eynrm2 — kind10: returns ||X||₂ for complex X (real result, complex input). */
#include <math.h>
typedef _Complex long double T;
typedef long double R;

R eynrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    if (n < 1 || incx < 1) return 0.0L;
    if (n == 1) {
        R re = __real__ x[0], im = __imag__ x[0];
        return hypotl(re, im);
    }

    /* Pass 1: max(|re|, |im|) → scale. */
    R scale = 0.0L;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        R re = fabsl(__real__ x[ix]), im = fabsl(__imag__ x[ix]);
        if (re > scale) scale = re;
        if (im > scale) scale = im;
        ix += incx;
    }
    if (scale == 0.0L) return 0.0L;
    if (!isfinite(scale)) return scale;

    /* Pass 2: (re/scale)² + (im/scale)². */
    R s = 0.0L;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        R re = __real__ x[ix] / scale, im = __imag__ x[ix] / scale;
        s += re * re + im * im;
        ix += incx;
    }
    return scale * sqrtl(s);
}
