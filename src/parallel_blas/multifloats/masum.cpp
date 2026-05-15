/* masum — multifloats real DD: returns Σ |X|. */
#include <cstddef>
#include <multifloats.h>
#include <multifloats/float64x2.h>

namespace mf = multifloats;
using T = mf::float64x2;

extern "C" T masum_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    T s{0.0, 0.0};
    if (n < 1 || incx < 1) return s;
    T s0{0.0, 0.0}, s1{0.0, 0.0};
    if (incx == 1) {
        int i = 0;
        for (; i + 1 < n; i += 2) {
            s0 = s0 + fabsdd(x[i]);
            s1 = s1 + fabsdd(x[i + 1]);
        }
        if (i < n) s0 = s0 + fabsdd(x[i]);
    } else {
        for (int i = 0, ix = 0; i < n; ++i, ix += incx)
            s0 = s0 + fabsdd(x[ix]);
    }
    return s0 + s1;
}
