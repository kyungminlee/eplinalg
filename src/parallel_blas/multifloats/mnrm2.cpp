/* mnrm2 — multifloats real DD: returns ||X||₂ via two-pass scaled. */
#include <cstddef>
#include <multifloats.h>
#include <multifloats/float64x2.h>

namespace mf = multifloats;
using T = mf::float64x2;

namespace {
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_gt(T a, T b) {
    return a.limbs[0] > b.limbs[0]
        || (a.limbs[0] == b.limbs[0] && a.limbs[1] > b.limbs[1]);
}
}

extern "C" T mnrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    T zero{0.0, 0.0};
    if (n < 1 || incx < 1) return zero;
    if (n == 1) return fabsdd(x[0]);

    /* Pass 1: scale = max|X|. */
    T scale = zero;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        T ax = fabsdd(x[ix]);
        if (dd_gt(ax, scale)) scale = ax;
        ix += incx;
    }
    if (dd_iszero(scale)) return zero;

    /* Pass 2: sum (X/scale)² then scale·√sum. */
    T s = zero;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        T t = x[ix] / scale;
        s = s + t * t;
        ix += incx;
    }
    return scale * sqrtdd(s);
}
