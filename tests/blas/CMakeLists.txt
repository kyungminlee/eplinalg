# Differential precision tests for migrated BLAS.
#
# Each test program runs the same operation through (a) the migrated
# extended-precision routine under test (linked from ${LIB_PAIR_PREFIX}blas)
# and (b) a REAL(KIND=16) reference (refblas_quad — see refblas/), then
# compares results in REAL(KIND=16) and writes a per-test JSON report.
#
# Reference library selection: -DBLAS_REF=netlib (default) builds a
# quad-precision reference from the vendored Netlib sources via
# gfortran's -freal-8-real-16; -DBLAS_REF=system links the system
# -lblas (precision capped at ~1e-16 in that mode).

# Quad precision is required: every test declares REAL(KIND=16).
if(NOT HAVE_REAL16)
    message(STATUS "tests/blas: skipping — HAVE_REAL16 is required")
    return()
endif()

set(BLAS_REF "netlib" CACHE STRING
    "Source for reference BLAS in tests/blas: netlib | system")
set_property(CACHE BLAS_REF PROPERTY STRINGS netlib system)

# ── Reference quad-precision BLAS ───────────────────────────────────
add_subdirectory(refblas)
if(NOT TARGET refblas_quad)
    message(STATUS "tests/blas: skipping — refblas_quad target not built")
    return()
endif()

# ── Common helper modules ───────────────────────────────────────────
add_library(blas_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    ../common/test_data.f90
    common/ref_quad_blas.f90
)
fortran_module_layout(blas_test_common)
target_link_libraries(blas_test_common PUBLIC refblas_quad)

# ── Per-target wrapper module ───────────────────────────────────────
# TARGET_NAME comes from target_config.cmake (kind10 / kind16 / multifloats)
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/blas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_blas_wrapper_src
    ${_wrapper_dir}/target_blas.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(blas_test_target STATIC ${_blas_wrapper_src})
fortran_module_layout(blas_test_target)
# Baseline (kind4/kind8): link the standard Netlib archive directly
# (s/c or d/z entry points), since there is no migrated archive to
# point at. Migrated targets (kind10/kind16/multifloats) link
# ${LIB_PAIR_PREFIX}blas as before.
if(BASELINE_BUILD)
    target_link_libraries(blas_test_target PUBLIC
        blas_test_common
        target_conv
        blas
    )
else()
    target_link_libraries(blas_test_target PUBLIC
        blas_test_common
        target_conv
        ${LIB_PAIR_PREFIX}blas
    )
endif()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(blas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
# Reports land in <build>/precision_reports/<test>.json, one JSON
# document per test program (multiple cases per file).
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

function(add_blas_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_link_libraries(${test_name} PRIVATE blas_test_target)
    add_test(NAME blas_${test_name}
             COMMAND ${test_name}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_lvl level1 level2 level3)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_lvl}/test_*.f90)
    foreach(_t ${_tests})
        add_blas_precision_test(${_t})
    endforeach()
endforeach()
