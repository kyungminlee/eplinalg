/* edot — kind10 real: returns Σ X·Y. */
typedef long double T;

T edot_(const int *n_, const T *x, const int *incx_,
        const T *y, const int *incy_)
{
    const int n = *n_, incx = *incx_, incy = *incy_;
    if (n <= 0) return 0.0L;
    T s = 0.0L;
    if (incx == 1 && incy == 1) {
        /* Two-accumulator unroll to mask x87 fadd latency. */
        T s0 = 0.0L, s1 = 0.0L;
        int i = 0;
        for (; i + 1 < n; i += 2) {
            s0 += x[i]     * y[i];
            s1 += x[i + 1] * y[i + 1];
        }
        s = s0 + s1;
        for (; i < n; ++i) s += x[i] * y[i];
    } else {
        int ix = (incx < 0) ? (-n + 1) * incx : 0;
        int iy = (incy < 0) ? (-n + 1) * incy : 0;
        for (int i = 0; i < n; ++i) { s += x[ix] * y[iy]; ix += incx; iy += incy; }
    }
    return s;
}
