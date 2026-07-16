/* multifloats_bridge.h -- C++ bridge for migrated BLACS / PBLAS sources.
 *
 * Exposes the upstream multifloats v0.4.0 names ``float64x2`` and
 * ``complex64x2`` (also the names emitted by the C migrator) at file
 * scope, and adds:
 *
 *   - mixed-type operators (``double op float64x2``)
 *   - MF_IS_ZERO / MF_IS_ONE macros for the C migrator's PBLAS rewrites
 *   - mf_abs / mf_cabs1 helpers (BLACS ABS / CABS1 macro replacements)
 *   - MPI datatype + op handle declarations (registered by
 *     multifloats_mpi_init, defined in multifloats_mpi.cpp)
 *
 * Usage: compile with  -include multifloats_bridge.h
 *
 * Requires: the multifloats v0.4.0 headers on the include path (either
 * installed, or from a FetchContent'd checkout). The upstream POD
 * ``multifloats::float64x2`` is the operator-rich C++ class we re-export
 * as global ``float64x2``; ``multifloats::complex64x2`` is the matching
 * POD struct with ``.re`` / ``.im`` members of type ``float64x2``.
 */
#ifndef MULTIFLOATS_BRIDGE_H
#define MULTIFLOATS_BRIDGE_H

#ifdef __cplusplus

/* PBLAS's ``pblas.h`` defines ``#define Int int`` at file scope, which
 * breaks multifloats's ``template <typename Int>`` parameter name.
 * Push/pop the macro around the upstream include so the template
 * parameter is parsed verbatim while leaving the rest of the TU's use
 * of ``Int`` alone.
 */
#pragma push_macro("Int")
#undef Int
#include <multifloats.h>
#pragma pop_macro("Int")

/* Skip the C++ MPI bindings (mpicxx.h). Without this guard, any migrated
 * source compiled as C++ that transitively pulls in <mpi.h> through this
 * bridge gets thousands of template declarations from mpicxx.h. Those
 * templates cannot live inside the ``extern "C" { … }`` wrap that the
 * c_migrator post-pass injects around .c bodies, so the link of
 * scalapack_c (whose REDIST sources include redist.h ->
 * multifloats_bridge.h -> mpi.h) fails. Defining MPICH_SKIP_MPICXX /
 * OMPI_SKIP_MPICXX before the include is the documented way to compile
 * MPI clients without the C++ bindings. */
#define MPICH_SKIP_MPICXX 1
#define OMPI_SKIP_MPICXX 1
#include <mpi.h>

/* ------------------------------------------------------------------ */
/* Bring the upstream types into the global namespace under the same  */
/* names the C migrator emits.                                        */
/* ------------------------------------------------------------------ */

using float64x2 = multifloats::float64x2;

/* The migrated BLACS sources reference complex members through the
 * BLACS DCOMPLEX convention (.r/.i), while upstream's POD spells them
 * .re/.im. Provide a layout-compatible local struct with anonymous
 * unions so both spellings resolve. The 4×double layout matches
 * ``multifloats::complex64x2``; reinterpret_cast between the two is
 * safe and used by the MPI op kernels. */
struct complex64x2 {
    union { float64x2 r; float64x2 re; };
    union { float64x2 i; float64x2 im; };
    complex64x2() : r{}, i{} {}
    complex64x2(float64x2 re_, float64x2 im_) : r(re_), i(im_) {}
};

/* PBLAS / BLACS C array form for double-double complex.
 * The migrator's ``c_type_aliases`` in pblas.yaml rewrites
 * ``cmplx`` / ``cmplx16`` to ``cmplx{RPU}`` (= ``cmplxM`` for the
 * current multifloats M-prefix; was ``cmplxT`` under the prior T/V
 * scheme and ``cmplxDD`` under the original DD/ZZ scheme).
 * All three names alias the same 2×float64x2 layout so existing
 * mfc_overrides referring to legacy spellings keep working. */
typedef float64x2 cmplxM[2];
typedef float64x2 cmplxT[2];
typedef float64x2 cmplxDD[2];

/* ------------------------------------------------------------------ */
/* MF_IS_* macros injected by c_migrator for C scalar comparisons.    */
/* In C++ they just use the native == operator on MultiFloat.         */
/* ------------------------------------------------------------------ */

#define MF_IS_ZERO(x) ((x) == float64x2{0.0})
#define MF_IS_ONE(x)  ((x) == float64x2{1.0})

/* ------------------------------------------------------------------ */
/* Mixed-type operators.                                              */
/*                                                                    */
/* The upstream class defines only same-type operators plus a         */
/* non-explicit constructor from ``double``, so ``x op 1.0`` works    */
/* via implicit conversion. The reversed form ``1.0 op x`` needs      */
/* these free functions — PBLAS uses idioms like ``ONE + ssq*...``    */
/* with ONE=1.0 (double).                                             */
/* ------------------------------------------------------------------ */

inline float64x2 operator+(double a, float64x2 b) { return float64x2{a} + b; }
inline float64x2 operator-(double a, float64x2 b) { return float64x2{a} - b; }
inline float64x2 operator*(double a, float64x2 b) { return float64x2{a} * b; }
inline float64x2 operator/(double a, float64x2 b) { return float64x2{a} / b; }
inline bool operator==(double a, float64x2 b) { return float64x2{a} == b; }
inline bool operator!=(double a, float64x2 b) { return float64x2{a} != b; }
inline bool operator<(double a, float64x2 b)  { return float64x2{a} < b; }
inline bool operator>(double a, float64x2 b)  { return float64x2{a} > b; }
inline bool operator<=(double a, float64x2 b) { return float64x2{a} <= b; }
inline bool operator>=(double a, float64x2 b) { return float64x2{a} >= b; }

/* ------------------------------------------------------------------ */
/* Integer truncation helper.                                          */
/* ------------------------------------------------------------------ */

inline int mf_to_int(float64x2 x) {
    return static_cast<int>(static_cast<double>(x));
}

/* ------------------------------------------------------------------ */
/* Absolute value (for BLACS ABS macro compatibility).                */
/* ------------------------------------------------------------------ */

inline float64x2 mf_abs(float64x2 x) { return x < float64x2{0.0} ? -x : x; }
inline float64x2 mf_cabs1(complex64x2 z) { return mf_abs(z.re) + mf_abs(z.im); }

/* ------------------------------------------------------------------ */
/* MPI handle externs + init (defined in multifloats_mpi.cpp).        */
/* ------------------------------------------------------------------ */

extern "C" {

extern MPI_Datatype MPI_FLOAT64X2;
extern MPI_Datatype MPI_COMPLEX64X2;

extern MPI_Op MPI_MM_SUM;
extern MPI_Op MPI_WW_SUM;
extern MPI_Op MPI_MM_AMX;
extern MPI_Op MPI_MM_AMN;
extern MPI_Op MPI_WW_AMX;
extern MPI_Op MPI_WW_AMN;

/* Deprecated pre-v0.14 names, kept one release cycle as aliases of the
 * MPI_MM_* / MPI_WW_* handles above (defined in multifloats_mpi.cpp). */
#if defined(__GNUC__) || defined(__clang__)
#  define MF_MPI_DEPRECATED(new_name) \
     __attribute__((deprecated("renamed to " #new_name " in v0.14")))
#else
#  define MF_MPI_DEPRECATED(new_name)
#endif
extern MPI_Op MPI_DD_SUM MF_MPI_DEPRECATED(MPI_MM_SUM);
extern MPI_Op MPI_DD_AMX MF_MPI_DEPRECATED(MPI_MM_AMX);
extern MPI_Op MPI_DD_AMN MF_MPI_DEPRECATED(MPI_MM_AMN);
extern MPI_Op MPI_ZZ_SUM MF_MPI_DEPRECATED(MPI_WW_SUM);
extern MPI_Op MPI_ZZ_AMX MF_MPI_DEPRECATED(MPI_WW_AMX);
extern MPI_Op MPI_ZZ_AMN MF_MPI_DEPRECATED(MPI_WW_AMN);

void multifloats_mpi_init(void);

} /* extern "C" */

#endif /* __cplusplus */
#endif /* MULTIFLOATS_BRIDGE_H */
