# Differential precision tests for migrated XBLAS.
#
# XBLAS routines compute the same mathematical operations as ordinary
# BLAS — the ``_x`` suffix only changes the internal accumulation
# precision used to reach the result. So our reference is the
# ``-freal-8-real-16``-promoted Netlib BLAS already present in
# ``refblas_quad`` (built by tests/blas/refblas), plus a few
# hand-coded quad-precision references for XBLAS-only routines
# (axpby, sum, waxpby, gemv2, gbmv2, symv2, hemv2, ge_sum_mv) — see
# ``common/ref_quad_xblas.f90``.
#
# Each test program runs the same operation through (a) the migrated
# extra-precise routine under test (linked from ${LIB_PREFIX}xblas)
# and (b) the quad reference, then compares results in REAL(KIND=16)
# and writes a per-test JSON report.

# Quad precision is required: every test declares REAL(KIND=16).
if(NOT HAVE_REAL16)
    message(STATUS "tests/xblas: skipping — HAVE_REAL16 is required")
    return()
endif()

# XBLAS tests reuse the quad-promoted Netlib BLAS reference built by
# tests/blas/refblas/. Without it there is no canonical numerical
# answer to compare against.
if(NOT TARGET refblas_quad)
    message(STATUS "tests/xblas: skipping — refblas_quad not available "
                   "(tests/blas must be enabled)")
    return()
endif()

# ── Per-target wrapper presence check ───────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/xblas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# tests/xblas/common reuses the same module names as tests/blas/common
# (prec_kinds, compare, prec_report, test_data); a separate .mod
# directory prevents collisions when both test trees build together.
set(_xb_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}/xblas_tests")
file(MAKE_DIRECTORY ${_xb_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(xblas_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    ../common/test_data.f90
    common/ref_quad_xblas.f90
)
set_target_properties(xblas_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_xb_moddir}"
)
target_include_directories(xblas_test_common PUBLIC
    $<BUILD_INTERFACE:${_xb_moddir}>
)
target_link_libraries(xblas_test_common PUBLIC refblas_quad)

# ── Per-target wrapper module ───────────────────────────────────────
add_fypp_source(_xblas_wrapper_src
    ${_wrapper_dir}/target_xblas.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(xblas_test_target STATIC ${_xblas_wrapper_src})
set_target_properties(xblas_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_xb_moddir}"
)
target_include_directories(xblas_test_target PUBLIC
    $<BUILD_INTERFACE:${_xb_moddir}>
)
target_link_libraries(xblas_test_target PUBLIC
    xblas_test_common
    target_conv
    ${LIB_PREFIX}xblas
    ${LIB_PREFIX}blas
)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(xblas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

function(add_xblas_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE xblas_test_target)
    add_test(NAME xblas_${test_name}
             COMMAND ${test_name}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_lvl level1 level2 level3)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_lvl}/test_*.f90)
    foreach(_t ${_tests})
        add_xblas_precision_test(${_t})
    endforeach()
endforeach()
