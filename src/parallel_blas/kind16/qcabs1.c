/* qcabs1 — kind16: |re(z)| + |im(z)|. */
#include <quadmath.h>
typedef __complex128 T;
typedef __float128 R;
R qcabs1_(const T *z) { return fabsq(__real__ *z) + fabsq(__imag__ *z); }
