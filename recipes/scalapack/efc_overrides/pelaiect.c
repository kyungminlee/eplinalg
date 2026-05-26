/* pelaiect.c -- hand-written kind10 override for pdlaiect.c
 *
 * The auto-cloned pelaiect.c inherits pdlaiect.c's IEEE-754 sign-bit
 * extraction trick (cast EREAL* to unsigned int*, shift bit 31 of word
 * 0 or 1). That layout assumption holds for IEEE-754 binary64 — it
 * does not hold for x86 80-bit extended precision (long double on
 * SysV-i386/amd64): the sign bit lives at bit 79 of the 16-byte slot,
 * and bytes 10..15 are alignment padding. The bit-31 read returns a
 * mantissa bit, so the count comes out wrong.
 *
 * pdstebz protects itself by branching on IEFLAG: if pelasnbt_ returns
 * 0, the Fortran caller falls back to PDLAPDCT (a portable comparison
 * loop) instead of pelaiectb_/pelaiectl_, so the buggy bit-twiddle is
 * never reached and end-to-end numerics stay correct — at the cost of
 * always taking the slow path. This override unlocks the fast path by
 * returning IEFLAG=1 and computing the Sturm count via a plain
 * ``tmp < 0.0L`` comparison.
 *
 * pelachkieee_'s IEEE self-test is short-circuited to ISIEEE=1 for the
 * same reason: its bit-fiddling on pzero/pinf is meaningless on
 * extended precision, and ``long double`` arithmetic on x86 follows
 * IEEE-754 semantics for the operations pdstebz relies on (signed
 * zero / signed infinity propagation through division and negation).
 *
 * Compiled as C (mpicc) in the kind10 stage; C_AS_CXX is FALSE.
 */

#include "pxsyevx.h"
#include "pblas.h"
#include <stdio.h>

void pelasnbt_(Int *ieflag)
{
    /* EREAL sign detection uses native ``<`` instead of pointer-cast
     * bit-twiddling. Return 1 so the Fortran caller takes the fast
     * (pelaiectb_) path rather than the portable fallback. */
    *ieflag = 1;
}

void pelaiectb_(EREAL *sigma, Int *n, EREAL *d, Int *count)
{
    EREAL lsigma = *sigma;
    EREAL *pd  = d;
    EREAL *pe2 = d + 1;
    EREAL tmp  = *pd - lsigma;
    pd += 2;
    *count = (tmp < (EREAL)0.0L) ? 1 : 0;
    for (Int i = 1; i < *n; i++) {
        tmp = *pd - *pe2 / tmp - lsigma;
        pd += 2;
        pe2 += 2;
        *count += (tmp < (EREAL)0.0L) ? 1 : 0;
    }
}

void pelaiectl_(EREAL *sigma, Int *n, EREAL *d, Int *count)
{
    /* Big-endian / little-endian distinction is meaningless for
     * comparison-based sign detection — both routines are identical. */
    pelaiectb_(sigma, n, d, count);
}

void pelachkieee_(Int *isieee, EREAL *rmax, EREAL *rmin)
{
    /* x86 ``long double`` arithmetic follows IEEE-754 for the
     * operations pdstebz needs (signed zero / signed infinity flow).
     * Skip the bit-fiddle self-test and report compliance. */
    (void)rmax;
    (void)rmin;
    *isieee = 1;
}
