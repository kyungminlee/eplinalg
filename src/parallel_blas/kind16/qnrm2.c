/* qnrm2 — kind16 real: returns ||X||₂. Two-pass scaled. */
#include <quadmath.h>
/* fabsq via __builtin_fabsf128 — single `pand` instead of a libquadmath function call. */
#undef fabsq
#define fabsq(x) __builtin_fabsf128(x)
typedef __float128 T;

T qnrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    if (n < 1 || incx < 1) return 0.0Q;
    if (n == 1) return fabsq(x[0]);
    T scale = 0.0Q;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        T ax = fabsq(x[ix]);
        if (ax > scale) scale = ax;
        ix += incx;
    }
    if (scale == 0.0Q) return 0.0Q;
    if (!__builtin_isfinite(scale)) return scale;
    T s = 0.0Q;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        T t = x[ix] / scale;
        s += t * t;
        ix += incx;
    }
    return scale * sqrtq(s);
}
