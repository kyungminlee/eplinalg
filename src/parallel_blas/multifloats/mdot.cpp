/* mdot — multifloats real DD: returns Σ X·Y.
 *
 * Returns float64x2 (16 bytes) via System V xmm0:xmm1 register pair —
 * matches gfortran's TYPE(real64x2) FUNCTION convention since multifloats
 * uses bind(c) types.
 */
#include <cstddef>
#include <multifloats.h>

namespace mf = multifloats;
using T = mf::float64x2;

extern "C" T mdot_(const int *n_,
                   const T *x, const int *incx_,
                   const T *y, const int *incy_)
{
    const int n = *n_, incx = *incx_, incy = *incy_;
    T s{0.0, 0.0};
    if (n <= 0) return s;
    if (incx == 1 && incy == 1) {
        T s0{0.0, 0.0}, s1{0.0, 0.0};
        int i = 0;
        for (; i + 1 < n; i += 2) {
            s0 = s0 + x[i]     * y[i];
            s1 = s1 + x[i + 1] * y[i + 1];
        }
        s = s0 + s1;
        for (; i < n; ++i) s = s + x[i] * y[i];
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        int iy = (incy < 0) ? (-n + 1) * incy : 0;
        for (int i = 0; i < n; ++i) { s = s + x[ix] * y[iy]; ix += incx; iy += incy; }
    }
    return s;
}
