/* wdotu — multifloats complex DD: returns Σ X·Y (unconjugated).
 *
 * Return type is complex64x2 (32 bytes) — System V ABI returns via
 * hidden sret pointer; gcc handles this transparently with return-by-value.
 */
#include <cstddef>
#include <multifloats.h>

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {
inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
}

extern "C" T wdotu_(const int *n_,
                    const T *x, const int *incx_,
                    const T *y, const int *incy_)
{
    const int n = *n_, incx = *incx_, incy = *incy_;
    T s{R{0.0, 0.0}, R{0.0, 0.0}};
    if (n <= 0) return s;
    if (incx == 1 && incy == 1) {
        for (int i = 0; i < n; ++i) s = cadd(s, cmul(x[i], y[i]));
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        int iy = (incy < 0) ? (-n + 1) * incy : 0;
        for (int i = 0; i < n; ++i) { s = cadd(s, cmul(x[ix], y[iy])); ix += incx; iy += incy; }
    }
    return s;
}
