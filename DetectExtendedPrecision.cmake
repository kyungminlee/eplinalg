# Detect compiler support for 80-bit and 128-bit extended real types,
# emitting ``add_compile_definitions(HAVE_REAL10)`` /
# ``add_compile_definitions(HAVE_REAL16)`` for each that probes
# successfully. The migrated ``la_constants_ep.F90`` /
# ``la_xisnan_ep.F90`` gate their qp/ep blocks on these macros so the
# build stays green on compilers that only support one of the two
# extended types (notably LLVM Flang 20, which has no working
# REAL(KIND=16)).
#
# Two CMake driver paths include this file:
#   - cmake/CMakeLists.txt           (shared, used by ``migrator stage``)
#   - the embedded template in       (per-library, used by ``migrator run``)
#     src/migrator/__main__.py       — copies this file alongside.
include_guard(GLOBAL)

include(CheckFortranSourceCompiles)

# 80-bit (KIND=10). Probe the literal ``REAL(KIND=10)`` on purpose,
# not ``selected_real_kind(18, 4931)``. On compilers that lack x87
# 80-bit reals (ifx on Linux/x86-64), the latter returns 16 and
# silently aliases the KIND=10 path onto KIND=16, which makes the
# EISNAN and QISNAN overloads in LA_XISNAN_EP collapse to the same
# argument type and triggers an ambiguous-generic-interface error.
# Asking for KIND=10 directly fails to compile on such compilers,
# which is the behavior we want.
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
# (the pattern used inside la_constants_ep to pick ``qp``), not the
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
