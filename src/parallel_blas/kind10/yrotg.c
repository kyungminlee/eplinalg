/* yrotg — kind10 complex Givens generator (ZROTG analog). */
#include <math.h>
typedef _Complex long double T;
typedef long double R;

void yrotg_(T *a_, const T *b_, R *c, T *s)
{
    const T a = *a_, b = *b_;
    R ar = __real__ a, ai = __imag__ a;
    R br = __real__ b, bi = __imag__ b;
    R aa = hypotl(ar, ai);
    R bb = hypotl(br, bi);
    if (aa == 0.0L) { *c = 0.0L; *s = 1.0L; *a_ = b; return; }
    if (bb == 0.0L) { *c = 1.0L; *s = 0.0L; return; }
    R t = hypotl(aa, bb);
    *c = aa / t;
    /* sign(a) = a / |a| */
    R sr = ar / aa, si = ai / aa;
    /* s = sign(a) · conj(b) / t */
    T sgn = sr + si * 1.0iL;
    T cb = br - bi * 1.0iL;
    *s = sgn * cb / t;
    *a_ = sgn * t;
}
