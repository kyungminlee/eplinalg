/* enrm2 — kind10 real: returns ||X||₂ = sqrt(Σ X·X).
 *
 * Simple two-pass scaled implementation: find max(|X|), then sum
 * (X/scale)². Avoids overflow on the inner sum. Doesn't use Blue's
 * algorithm but handles the dynamic-range case correctly.
 */
#include <math.h>
typedef long double T;

T enrm2_(const int *n_, const T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    if (n < 1 || incx < 1) return 0.0L;
    if (n == 1) return fabsl(x[0]);

    /* Pass 1: find scale = max|X|. */
    T scale = 0.0L;
    int ix = 0;
    for (int i = 0; i < n; ++i) {
        T ax = fabsl(x[ix]);
        if (ax > scale) scale = ax;
        ix += incx;
    }
    if (scale == 0.0L) return 0.0L;
    if (!isfinite(scale)) return scale;  /* propagate Inf/NaN */

    /* Pass 2: sum (X/scale)². */
    T s = 0.0L;
    ix = 0;
    for (int i = 0; i < n; ++i) {
        T t = x[ix] / scale;
        s += t * t;
        ix += incx;
    }
    return scale * sqrtl(s);
}
