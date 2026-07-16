/*
 * Target-aware real- and complex-arithmetic shim for the C-side MUMPS
 * tests.
 *
 * The C tests (test_*_c_basic.c, test_*_c_sym.c) need a real type wide
 * enough to interoperate with the migrated bridge's typedef of
 * DMUMPS_REAL.  Each precision target ships its own
 * tests/mumps/target_<TARGET_NAME>/c/include/mumps_c_types.h that
 * defines DMUMPS_REAL appropriately:
 *
 *   kind16       — __float128
 *   kind10       — long double
 *   multifloats  — POD ``mumps_float64x2`` (no scalar arithmetic in C)
 *
 * On kind16 / kind10 the test side uses the same type as the bridge —
 * ``test_real`` is the precision-wide scalar and ``test_complex`` is
 * the precision-wide ``{r, i}`` struct.  TR_* macros bind to the
 * matching literal suffix, math helpers and snprintf wrapper.
 *
 * On multifloats the bridge type is a POD struct with no plain-C
 * arithmetic.  The test side stays in plain ``double`` (host
 * precision) and copies through narrow ↔ wide buffers at the bridge
 * boundary; ``tr_widen`` / ``tr_narrow`` (real) and ``tc_widen`` /
 * ``tc_narrow`` (complex) translate between ``test_real`` /
 * ``test_complex`` and ``mumps_float64x2`` / ``mumps_complex64x2``.
 * The 4×4 hand-built test problems don't strain double precision; the
 * point of the C tests is to verify the bridge wires correctly to the
 * multifloats solver, not to test multifloats arithmetic in C.
 *
 * Selected via the ``TEST_TARGET_<NAME>`` define passed by the cmake
 * glue (see ``add_mumps_c_test`` in tests/mumps/CMakeLists.txt).
 */

#ifndef TEST_REAL_COMPAT_H
#define TEST_REAL_COMPAT_H

#include <stdio.h>

/* Pull in the per-target real/complex types so test_complex can be
 * typedef'd to the target's own extended complex type (mumps_float128_complex
 * on kind16, mumps_long_double_complex on kind10 — both defined by the
 * migrated bridge header this test #includes), and mumps_float64x2 /
 * mumps_complex64x2 are visible on multifloats. Using the bridge's exact
 * complex type (rather than the double-precision mumps_double_complex) keeps
 * the complex C test at the target's true precision AND lets id.a / id.rhs
 * take test_complex arrays directly with no type mismatch. */
#include TARGET_REAL_HEADER
#include TARGET_COMPLEX_HEADER

#if defined(TEST_TARGET_KIND16)
#  include <quadmath.h>
   typedef __float128 test_real;
   typedef mumps_float128_complex test_complex;
#  define TR_LIT(x)      x##q
#  define TR_FABS(x)     fabsq(x)
#  define TR_SQRT(x)     sqrtq(x)
#  define TR_EPS         FLT128_EPSILON
#  define TR_MIN         FLT128_MIN
   static inline int test_real_snprintf(char *buf, size_t n, test_real x) {
       return quadmath_snprintf(buf, n, "%.6Qe", x);
   }
#elif defined(TEST_TARGET_KIND10)
#  include <float.h>
#  include <math.h>
   typedef long double test_real;
   typedef mumps_long_double_complex test_complex;
#  define TR_LIT(x)      x##L
#  define TR_FABS(x)     fabsl(x)
#  define TR_SQRT(x)     sqrtl(x)
#  define TR_EPS         LDBL_EPSILON
#  define TR_MIN         LDBL_MIN
   static inline int test_real_snprintf(char *buf, size_t n, test_real x) {
       return snprintf(buf, n, "%.6Le", x);
   }
#elif defined(TEST_TARGET_MULTIFLOATS)
#  include <float.h>
#  include <math.h>
   typedef double test_real;
   typedef struct { double r, i; } test_complex;
#  define TR_LIT(x)      (x)
#  define TR_FABS(x)     fabs(x)
#  define TR_SQRT(x)     sqrt(x)
#  define TR_EPS         DBL_EPSILON
#  define TR_MIN         DBL_MIN
   static inline int test_real_snprintf(char *buf, size_t n, test_real x) {
       return snprintf(buf, n, "%.6e", x);
   }
   /* Bridge-boundary widening: test_real (double) ↔ mumps_float64x2.
    * Lower limb is zero — the host data is plain double, so the
    * upper limb captures the value exactly and the lower limb is
    * unused.  Narrow sums the limbs (lower is zero in our path; the
    * sum keeps the test driver's read of id.rhs faithful even if
    * the solver populated the lower limb). */
   static inline mumps_float64x2 tr_widen(test_real x) {
       mumps_float64x2 r;
       r.limbs[0] = (double)x;
       r.limbs[1] = 0.0;
       return r;
   }
   static inline test_real tr_narrow(mumps_float64x2 x) {
       return (test_real)(x.limbs[0] + x.limbs[1]);
   }
   static inline mumps_complex64x2 tc_widen(test_complex z) {
       mumps_complex64x2 r;
       r.r.limbs[0] = z.r; r.r.limbs[1] = 0.0;
       r.i.limbs[0] = z.i; r.i.limbs[1] = 0.0;
       return r;
   }
   static inline test_complex tc_narrow(mumps_complex64x2 z) {
       test_complex r;
       r.r = z.r.limbs[0] + z.r.limbs[1];
       r.i = z.i.limbs[0] + z.i.limbs[1];
       return r;
   }
#else
#  error "test_real_compat.h: define TEST_TARGET_KIND16, TEST_TARGET_KIND10, or TEST_TARGET_MULTIFLOATS"
#endif

/* Per-target one-time MPI setup, called by each C test immediately after
 * MPI_Init and before the first MUMPS_C. This mirrors the Fortran
 * target_qmumps wrapper's lazy init:
 *   kind16      — registers the custom quad reduce ops (the MPI_QQ and
 *                 MPI_XX families); Intel MPI has no 16-byte-real reduce
 *                 kernel, so MUMPS's
 *                 first Allreduce SIGSEGVs / reports Invalid MPI_Op at
 *                 np >= 1 without them.
 *   multifloats — registers the custom datatypes + reduce ops.
 *   kind10      — none needed (long-double reductions use builtin ops).
 * Routing through this hook keeps the C tests target-agnostic (they call
 * test_target_mpi_init() regardless of precision). */
#if defined(TEST_TARGET_KIND16)
   extern void quad_mpi_init(void);
   static inline void test_target_mpi_init(void) { quad_mpi_init(); }
#elif defined(TEST_TARGET_MULTIFLOATS)
   extern void multifloats_mpi_init(void);
   static inline void test_target_mpi_init(void) { multifloats_mpi_init(); }
#else
   static inline void test_target_mpi_init(void) { }
#endif

#ifndef TEST_TARGET_NAME
#  error "test_real_compat.h: define TEST_TARGET_NAME (e.g. \"kind16\")"
#endif

#endif /* TEST_REAL_COMPAT_H */
