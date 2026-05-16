/* yscal — kind10 complex: X := α·X (α and X both complex).
 *
 * Unit-stride loop uses pointer-compare termination so the inner loop
 * matches gfortran's codegen for the migrated reference. The original
 * `for (int i = 0; i < n; ++i) x[i] *= alpha` kept a separate i-counter
 * that added 2 instructions per iter to an already-deep x87 complex-mul
 * body — ~13% slowdown.
 */
typedef _Complex long double T;

void yscal_(const int *n_, const T *alpha_, T *x, const int *incx_)
{
    const int n = *n_, incx = *incx_;
    const T alpha = *alpha_;
    if (n <= 0) return;
    if (alpha == (T)1.0L) return;
    if (incx == 1) {
        T *end = x + n;
        for (T *p = x; p < end; ++p) *p *= alpha;
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        for (int i = 0; i < n; ++i) { x[ix] *= alpha; ix += incx; }
    }
}
