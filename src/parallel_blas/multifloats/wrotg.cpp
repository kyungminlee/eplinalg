/* wrotg — multifloats: complex Givens generator (ZROTG analog).
 *   a, b complex; c real; s complex.
 *   c = |a| / sqrt(|a|² + |b|²)
 *   s = (a/|a|) · conj(b) / sqrt(|a|² + |b|²)
 *   r = (a/|a|) · sqrt(|a|² + |b|²);  a := r
 */
#include <multifloats.h>
#include <multifloats/float64x2.h>

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {
inline bool dd_eq0(R x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool cdd_eq0(T const &z) { return dd_eq0(z.re) && dd_eq0(z.im); }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }
inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cdiv_real(T const &a, R const &r) { return T{ a.re / r, a.im / r }; }
}

extern "C" void wrotg_(T *a_, const T *b_, R *c, T *s)
{
    const T a = *a_, b = *b_;
    const T czero{R{0.0, 0.0}, R{0.0, 0.0}};
    const T cone {R{1.0, 0.0}, R{0.0, 0.0}};
    if (cdd_eq0(a)) {
        *c = R{0.0, 0.0};
        *s = cone;
        *a_ = b;
        return;
    }
    R aa = cabsdd(a);          /* |a| (overflow-safe hypot) */
    if (cdd_eq0(b)) {
        *c = R{1.0, 0.0};
        *s = czero;
        return;
    }
    R bb = cabsdd(b);
    /* t = sqrt(|a|² + |b|²) — use hypot-like via cabs of a 2-component */
    R t = sqrtdd(aa * aa + bb * bb);
    *c = aa / t;
    /* sign(a) = a / |a| */
    T sgn = cdiv_real(a, aa);
    *s = cdiv_real(cmul(sgn, cconj(b)), t);
    *a_ = T{ sgn.re * t, sgn.im * t };
}
