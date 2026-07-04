# Differential precision tests for migrated LAPACK.
#
# Mirrors tests/blas/CMakeLists.txt. Each test program runs the same
# operation through (a) the migrated extended-precision routine under
# test (linked from ${LIB_PREFIX}lapack) and (b) a REAL(KIND=16)
# reference (reflapack_quad — see reflapack/, which in turn links
# refblas_quad from the BLAS side), then compares results in
# REAL(KIND=16) and writes a per-test JSON report.
#
# Reference library selection:
#   -DLAPACK_REF=netlib (default) builds a quad-precision reference
#   from the vendored Netlib LAPACK sources via gfortran's
#   -freal-8-real-16; -DLAPACK_REF=system links the system -llapack
#   (precision capped at ~1e-16 in that mode).

# Quad precision is required: every test declares REAL(KIND=16).
if(NOT HAVE_REAL16)
    message(STATUS "tests/lapack: skipping — HAVE_REAL16 is required")
    return()
endif()

# LAPACK tests depend on BLAS infrastructure (refblas_quad for the
# reference, blas_test_target's linkage is independent but the quad
# BLAS reference library is defined by tests/blas/refblas).
if(NOT TARGET refblas_quad)
    message(STATUS "tests/lapack: skipping — refblas_quad not available "
                   "(tests/blas must be enabled)")
    return()
endif()

set(LAPACK_REF "netlib" CACHE STRING
    "Source for reference LAPACK in tests/lapack: netlib | system")
set_property(CACHE LAPACK_REF PROPERTY STRINGS netlib system)

# ── Reference quad-precision LAPACK ─────────────────────────────────
add_subdirectory(reflapack)
if(NOT TARGET reflapack_quad)
    message(STATUS "tests/lapack: skipping — reflapack_quad target not built")
    return()
endif()

# ── Per-target wrapper presence check ───────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/lapack: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# tests/lapack/common reuses the same module names as tests/blas/common
# (prec_kinds, compare, prec_report, test_data). Placing the LAPACK
# variants in their own subdirectory under fmod prevents .mod-file
# collisions when both test trees build together.
set(_lp_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/lapack_tests")
file(MAKE_DIRECTORY ${_lp_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(lapack_test_common STATIC
    common/prec_kinds.f90
    common/compare.f90
    common/prec_report.f90
    common/test_data.f90
    common/ref_quad_lapack.f90
)
set_target_properties(lapack_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_lp_moddir}"
)
target_include_directories(lapack_test_common PUBLIC
    $<BUILD_INTERFACE:${_lp_moddir}>
)
target_link_libraries(lapack_test_common PUBLIC reflapack_quad)

# ── Per-target wrapper module ───────────────────────────────────────
add_fypp_source(_lapack_wrapper_src
    ${_wrapper_dir}/target_lapack.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(lapack_test_target STATIC ${_lapack_wrapper_src})
set_target_properties(lapack_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_lp_moddir}"
)
target_include_directories(lapack_test_target PUBLIC
    $<BUILD_INTERFACE:${_lp_moddir}>
)
if(BASELINE_BUILD)
    # kind4/kind8: link the standard Netlib lapack/blas archives. There
    # is no ${LIB_PREFIX}lapack target — LIB_PREFIX is empty under
    # baseline staging, so the migrated archives never get built.
    target_link_libraries(lapack_test_target PUBLIC
        lapack_test_common
        target_conv
        lapack
        blas
    )
else()
    target_link_libraries(lapack_test_target PUBLIC
        lapack_test_common
        target_conv
        ${LIB_PREFIX}lapack
        ${LIB_PREFIX}blas
    )
endif()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(lapack_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
# Reports land in <build>/precision_reports/<test>.<target>.json,
# co-mingling with the BLAS reports for the Python aggregator.
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

function(add_lapack_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE lapack_test_target)
    add_test(NAME lapack_${test_name}
             COMMAND ${test_name}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_group linear_solve factorization eigenvalue svd auxiliary)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_group}/test_*.f90)
    if(BASELINE_BUILD)
        # XBLAS-dependent iterative-refinement tests (xx-suffix and
        # rfsx-suffix routines: dgbsvxx, dgerfsx, zporfsx, …) need a
        # native-precision XBLAS archive that this build doesn't ship.
        # The migrated path links qxblas + a quad bridge in
        # reflapack_quad; for kind4/kind8 there is neither, so the
        # corresponding tests would fail to link. Skip them.
        list(FILTER _tests EXCLUDE REGEX "/test_[a-z]+(rfsx|svxx)\\.f90$")
    endif()
    foreach(_t ${_tests})
        add_lapack_precision_test(${_t})
    endforeach()
endforeach()
