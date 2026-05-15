/* easum — kind10 real: returns Σ |X|. */
#include <math.h>
typedef long double T;

T easum_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    if (n < 1 || incx < 1) return 0.0L;
    T s0 = 0.0L, s1 = 0.0L;
    if (incx == 1) {
        int i = 0;
        for (; i + 1 < n; i += 2) { s0 += fabsl(x[i]); s1 += fabsl(x[i + 1]); }
        if (i < n) s0 += fabsl(x[i]);
    } else {
        for (int i = 0, ix = 0; i < n; ++i, ix += incx) s0 += fabsl(x[ix]);
    }
    return s0 + s1;
}
