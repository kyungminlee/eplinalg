/* qxnrm2 — kind16: ||X||₂ for complex X (real result). Two-pass scaled. */
#include <quadmath.h>
typedef __complex128 T;
typedef __float128 R;

R qxnrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    if (n < 1 || incx < 1) return 0.0Q;
    if (n == 1) return hypotq(__real__ x[0], __imag__ x[0]);
    R scale = 0.0Q;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        R re = fabsq(__real__ x[ix]), im = fabsq(__imag__ x[ix]);
        if (re > scale) scale = re;
        if (im > scale) scale = im;
        ix += incx;
    }
    if (scale == 0.0Q) return 0.0Q;
    if (!finiteq(scale)) return scale;
    R s = 0.0Q;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        R re = __real__ x[ix] / scale, im = __imag__ x[ix] / scale;
        s += re * re + im * im;
        ix += incx;
    }
    return scale * sqrtq(s);
}
