/* xrotg — kind16 complex Givens generator. */
#include <quadmath.h>
typedef __complex128 T;
typedef __float128 R;

void xrotg_(T *a_, const T *b_, R *c, T *s)
{
    const T a = *a_, b = *b_;
    R ar = __real__ a, ai = __imag__ a;
    R br = __real__ b, bi = __imag__ b;
    R aa = hypotq(ar, ai);
    R bb = hypotq(br, bi);
    if (aa == 0.0Q) { *c = 0.0Q; *s = b; *a_ = b; return; }
    if (bb == 0.0Q) { *c = 1.0Q; __real__ *s = 0; __imag__ *s = 0; return; }
    R t = hypotq(aa, bb);
    *c = aa / t;
    R sr = ar / aa, si = ai / aa;
    T sgn; __real__ sgn = sr; __imag__ sgn = si;
    T cb; __real__ cb = br; __imag__ cb = -bi;
    *s = sgn * cb / t;
    *a_ = sgn * t;
}
