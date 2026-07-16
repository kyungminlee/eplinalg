/* XBLAS extra-precise BLAS stubs for kind4 / kind8 baseline builds.
 *
 * Standard LAPACK ships iterative-refinement helpers (dla_*_extended.f
 * and zla_*_extended.f) that call into Netlib XBLAS via a fixed set of
 * `blas_*_x_` entry points. The migrated test build resolves those
 * via the ${LIB_PREFIX}xblas archive plus a quad-precision bridge in
 * tests/lapack/reflapack/refxblas_quad_bridge.f90; the kind4 / kind8
 * baseline build has no XBLAS at all.
 *
 * The dla_*_extended routines themselves are pulled into the test
 * executable transitively through target_lapack's public wrappers
 * (dgerfsx, dgesvxx, dgbrfsx, …) even when those baseline tests are
 * excluded — ld's lazy archive extraction still demands link-time
 * resolution of every undefined reference inside the extracted .o
 * files. We satisfy those references with no-op stubs that abort if
 * actually called. Baseline tests never invoke an *xx / *rfsx
 * routine (those tests are filtered out in tests/lapack/CMakeLists.txt
 * under BASELINE_BUILD), so the stubs stay quiescent at runtime.
 */

#include <stdio.h>
#include <stdlib.h>

static void xblas_baseline_unreachable(const char *name) {
    fprintf(stderr,
            "%s: XBLAS routine called from a kind4/kind8 baseline build.\n"
            "  Iterative-refinement (xx-suffix / rfsx-suffix) routines need\n"
            "  a native-precision XBLAS archive that this build does not\n"
            "  ship. Those tests are filtered out of the baseline ctest\n"
            "  set; reaching this stub means a non-baseline path slipped\n"
            "  through.\n",
            name);
    abort();
}

#define STUB(name)                                                        \
    void name(void) { xblas_baseline_unreachable(#name); }

STUB(blas_dgbmv_x_)
STUB(blas_dgbmv2_x_)
STUB(blas_dgemv_x_)
STUB(blas_dgemv2_x_)
STUB(blas_dsymv_x_)
STUB(blas_dsymv2_x_)
STUB(blas_zgbmv_x_)
STUB(blas_zgbmv2_x_)
STUB(blas_zgemv_x_)
STUB(blas_zgemv2_x_)
STUB(blas_zhemv_x_)
STUB(blas_zhemv2_x_)
STUB(blas_zsymv_x_)
STUB(blas_zsymv2_x_)
