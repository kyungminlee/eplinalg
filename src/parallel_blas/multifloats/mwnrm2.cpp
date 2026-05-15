/* mwnrm2 — multifloats: ||X||₂ for complex DD X, returns real DD. */
#include <cstddef>
#include <multifloats.h>
#include <multifloats/float64x2.h>

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {
inline bool dd_iszero(R x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_gt(R a, R b) {
    return a.limbs[0] > b.limbs[0]
        || (a.limbs[0] == b.limbs[0] && a.limbs[1] > b.limbs[1]);
}
}

extern "C" R mwnrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    R zero{0.0, 0.0};
    if (n < 1 || incx < 1) return zero;

    /* Pass 1: scale = max(|re|, |im|). */
    R scale = zero;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        R r = fabsdd(x[ix].re), m = fabsdd(x[ix].im);
        if (dd_gt(r, scale)) scale = r;
        if (dd_gt(m, scale)) scale = m;
        ix += incx;
    }
    if (dd_iszero(scale)) return zero;

    /* Pass 2: sum (re/scale)² + (im/scale)². */
    R s = zero;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        R re = x[ix].re / scale, im = x[ix].im / scale;
        s = s + re * re + im * im;
        ix += incx;
    }
    return scale * sqrtdd(s);
}
