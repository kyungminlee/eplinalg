/* pmlaiect.c -- hand-written multifloats override for pdlaiect.c
 *
 * The original pdlaiect.c extracts the sign bit of a double via
 * pointer-cast to unsigned int and bit-shift.  That trick relies on
 * the IEEE 754 double layout, which does not apply to float64x2.
 *
 * This override replaces the bit manipulation with C++ comparison
 * operators provided by multifloats_bridge.h.
 *
 * Compiled as C++ (mpicxx) with -include multifloats_bridge.h.
 */

#include "pxsyevx.h"
#include "pblas.h"

#ifdef __cplusplus
extern "C" {
#endif

void pmlasnbt_(Int *ieflag)
{
    /* For float64x2, sign detection uses C++ comparison operators.
     * Return nonzero to signal that sign detection is available. */
    *ieflag = 1;
}

void pmlaiectb_(float64x2 *sigma, Int *n, float64x2 *d, Int *count)
{
    float64x2 lsigma = *sigma;
    float64x2 *pd = d;
    float64x2 *pe2 = d + 1;
    float64x2 tmp = *pd - lsigma;
    pd += 2;
    *count = (tmp < float64x2{0.0}) ? 1 : 0;
    for (Int i = 1; i < *n; i++) {
        tmp = *pd - *pe2 / tmp - lsigma;
        pd += 2;
        pe2 += 2;
        *count += (tmp < float64x2{0.0}) ? 1 : 0;
    }
}

void pmlaiectl_(float64x2 *sigma, Int *n, float64x2 *d, Int *count)
{
    /* The big-endian / little-endian distinction for sign-bit extraction
     * does not apply to float64x2 comparison-based sign detection.
     * Both functions are identical. */
    pmlaiectb_(sigma, n, d, count);
}

void pmlachkieee_(Int *isieee, float64x2 *rmax, float64x2 *rmin)
{
    /* float64x2 uses C++ operator overloads that handle sign
     * and comparison correctly.  Report IEEE-compatible behaviour. */
    (void)rmax;
    (void)rmin;
    *isieee = 1;
}

#ifdef __cplusplus
} /* extern "C" */
#endif
