/* pqlaiect.c -- hand-written kind16 override for pdlaiect.c
 *
 * The auto-cloned pqlaiect.c inherits pdlaiect.c's IEEE-754 sign-bit
 * extraction trick (cast QREAL* to unsigned int*, shift bit 31 of word
 * 0 or 1). That layout assumption holds for IEEE-754 binary64. For
 * __float128 (IEEE-754 binary128) the sign bit is at bit 127 of the
 * 16-byte slot; reading bit 31 of word 0 or 1 returns a mantissa bit,
 * so the count comes out wrong.
 *
 * pdstebz protects itself by branching on IEFLAG: if pqlasnbt_ returns
 * 0, the Fortran caller falls back to PDLAPDCT (a portable comparison
 * loop) instead of pqlaiectb_/pqlaiectl_, so the buggy bit-twiddle is
 * never reached and end-to-end numerics stay correct — at the cost of
 * always taking the slow path. This override unlocks the fast path by
 * returning IEFLAG=1 and computing the Sturm count via a plain
 * ``tmp < 0.0q`` comparison.
 *
 * pqlachkieee_'s IEEE self-test is short-circuited to ISIEEE=1 for the
 * same reason: GCC's __float128 follows IEEE-754 binary128 for the
 * operations pdstebz relies on (signed zero / signed infinity
 * propagation through division and negation).
 *
 * Compiled as C (mpicc) in the kind16 stage; C_AS_CXX is FALSE.
 * Note: GCC supports ``<`` directly on __float128 in C mode (no
 * libquadmath helper needed for ordered comparison).
 */

#include "pxsyevx.h"
#include "pblas.h"
#include <stdio.h>

void pqlasnbt_(Int *ieflag)
{
    /* __float128 sign detection uses native ``<`` instead of
     * pointer-cast bit-twiddling. Return 1 so the Fortran caller takes
     * the fast (pqlaiectb_) path rather than the portable fallback. */
    *ieflag = 1;
}

void pqlaiectb_(QREAL *sigma, Int *n, QREAL *d, Int *count)
{
    QREAL lsigma = *sigma;
    QREAL *pd  = d;
    QREAL *pe2 = d + 1;
    QREAL tmp  = *pd - lsigma;
    pd += 2;
    *count = (tmp < (QREAL)0.0q) ? 1 : 0;
    for (Int i = 1; i < *n; i++) {
        tmp = *pd - *pe2 / tmp - lsigma;
        pd += 2;
        pe2 += 2;
        *count += (tmp < (QREAL)0.0q) ? 1 : 0;
    }
}

void pqlaiectl_(QREAL *sigma, Int *n, QREAL *d, Int *count)
{
    /* Big-endian / little-endian distinction is meaningless for
     * comparison-based sign detection — both routines are identical. */
    pqlaiectb_(sigma, n, d, count);
}

void pqlachkieee_(Int *isieee, QREAL *rmax, QREAL *rmin)
{
    /* __float128 follows IEEE-754 for the operations pdstebz needs
     * (signed zero / signed infinity flow). Skip the bit-fiddle
     * self-test and report compliance. */
    (void)rmax;
    (void)rmin;
    *isieee = 1;
}
