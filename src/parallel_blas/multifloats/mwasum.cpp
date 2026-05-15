/* mwasum — multifloats: Σ (|re(X)| + |im(X)|) for complex DD X. */
#include <cstddef>
#include <multifloats.h>
#include <multifloats/float64x2.h>

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

extern "C" R mwasum_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    R s{0.0, 0.0};
    if (n < 1 || incx < 1) return s;
    if (incx == 1) {
        for (int i = 0; i < n; ++i) s = s + fabsdd(x[i].re) + fabsdd(x[i].im);
    } else {
        for (int i = 0, ix = 0; i < n; ++i, ix += incx)
            s = s + fabsdd(x[ix].re) + fabsdd(x[ix].im);
    }
    return s;
}
