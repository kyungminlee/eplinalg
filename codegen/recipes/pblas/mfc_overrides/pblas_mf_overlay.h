/* pblas_mf_overlay.h -- multifloats macro overrides for PBLAS files
 * that do NOT use the algorithm-selection cost model.
 *
 * PBLAS C entry points share the macros ZERO, ONE, ABS defined in
 * PBtools.h. Their definitions assume primitive double scalars and
 * therefore stop working once the migrated file's data types become
 * the float64x2 struct (where ``=`` and ``< 0`` and unary ``-``
 * are not defined).
 *
 * Files that DO use the cost model (Level-3 BLAS like pdgemm, pdsymm)
 * declare local doubles ABest, ACest, BCest, tmp1..tmp4 and use
 * those macros at primitive-double type, so the c_migrator's
 * protect-double pass keeps them. Those files are NOT given this
 * overlay.
 *
 * Files that do NOT use the cost model (vector reductions: amax,
 * asum, dot, nrm2, her, her2, herk) need the macros to operate on
 * float64x2 directly. Including this overlay rewrites the macros
 * to use libmfc primitives.
 */
#ifndef PBLAS_MF_OVERLAY_H
#define PBLAS_MF_OVERLAY_H 1

#include "multifloats_bridge.h"

#ifdef ZERO
#  undef ZERO
#endif
#ifdef ONE
#  undef ONE
#endif
#ifdef ABS
#  undef ABS
#endif

#define ZERO   ((float64x2){{0.0, 0.0}})
#define ONE    ((float64x2){{1.0, 0.0}})
#define ABS(x) (mf_abs(x))

#endif /* PBLAS_MF_OVERLAY_H */
