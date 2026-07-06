# Detect compiler support for 80-bit and 128-bit extended real types,
# setting the CMake cache variables ``HAVE_REAL10`` / ``HAVE_REAL16``
# (and emitting the matching ``add_compile_definitions``) for each that
# probes successfully. The cache variables gate which extended-precision
# targets and tests are enabled (see the tests/*/CMakeLists.txt guards),
# so the build stays green on compilers that only support one of the two
# extended types (notably LLVM Flang 20, which has no working
# REAL(KIND=16)).
#
# The migrated per-target modules do not #ifdef on these macros: each
# lives in its precision-specific archive (la_constants_ey/la_xisnan_ey
# in libelapack, la_constants_qx/la_xisnan_qx in libqlapack) and is
# compiled only when that target is built — which by definition has the
# precision. The compile definitions are kept for any preprocessed
# source that wants a capability macro.
#
# Two CMake driver paths include this file:
#   - cmake/CMakeLists.txt           (shared, used by ``migrator stage``)
#   - the generated per-library CMakeLists from ``migrator run``
#     (``_generate_cmake`` in src/migrator/cmake_gen.py — copies this
#     file alongside).
include_guard(GLOBAL)

include(CheckFortranSourceCompiles)

# 80-bit (KIND=10). Probe the literal ``REAL(KIND=10)`` on purpose,
# not ``selected_real_kind(18, 4931)``. On compilers that lack x87
# 80-bit reals (ifx on Linux/x86-64), the latter returns 16 and
# silently aliases the KIND=10 path onto KIND=16. Asking for KIND=10
# directly fails to compile on such compilers, so HAVE_REAL10 stays
# undefined and LA_XISNAN_EY's body is #ifdef'd out entirely — which
# is the behavior we want. (The EY/QX split puts EISNAN and QISNAN in
# separate modules, so a KIND alias cannot produce an ambiguous generic
# interface between them.)
check_fortran_source_compiles("
    program test_real10
    real(kind=10) :: x
    x = 1.0_10
    end program
" HAVE_REAL10 SRC_EXT F90)
if(HAVE_REAL10)
    add_compile_definitions(HAVE_REAL10)
endif()

# 128-bit (IEEE binary128). Probe ``selected_real_kind(33, 4931)``
# (the pattern used inside la_constants_qx to pick ``qp``), not the
# literal ``REAL(KIND=16)``. LLVM Flang 20 accepts the KIND=16
# literal but then returns -1 for ``selected_real_kind(33)``, so the
# module ends up with ``qp = -1`` and every subsequent ``real(qp)``
# declaration fails. Probing the real selector catches this case.
check_fortran_source_compiles("
    program test_real16
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp) :: x
    x = 1.0_qp
    end program
" HAVE_REAL16 SRC_EXT F90)
if(HAVE_REAL16)
    add_compile_definitions(HAVE_REAL16)
endif()
