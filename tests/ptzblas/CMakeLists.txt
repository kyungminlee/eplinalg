# Differential precision tests for migrated PTZBLAS.
#
# PTZBLAS is the serial trapezoidal/auxiliary BLAS-like helper layer
# that PBLAS sits on top of — pure Fortran, no MPI, no BLACS. Each
# test program builds a small problem at REAL(KIND=16), calls the
# migrated single-letter-prefix entry point through a per-target
# wrapper, and compares against either refblas_quad (for the BLAS
# analogues — vasum/vvdot/dzvasum) or an inline quad-precision
# reference implementation in ptzblas_ref_quad_aux.

if(NOT HAVE_REAL16)
    message(STATUS "tests/ptzblas: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/ptzblas: skipping — refblas_quad target not built "
                   "(tests/blas must be staged alongside tests/ptzblas)")
    return()
endif()

if(NOT TARGET ${LIB_PAIR_PREFIX}ptzblas)
    message(STATUS "tests/ptzblas: skipping — ${LIB_PAIR_PREFIX}ptzblas not built")
    return()
endif()

# ── Common helper modules ───────────────────────────────────────────
add_library(ptzblas_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    ../common/test_data.f90
    common/ref_quad_blas.f90
    common/ref_quad_aux.f90
)
fortran_test_module_layout(ptzblas_test_common)
target_link_libraries(ptzblas_test_common PUBLIC refblas_quad)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/ptzblas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_ptzblas_wrapper_src
    ${_wrapper_dir}/target_ptzblas.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(ptzblas_test_target STATIC ${_ptzblas_wrapper_src})
fortran_test_module_layout(ptzblas_test_target)
target_link_libraries(ptzblas_test_target PUBLIC
    ptzblas_test_common
    target_conv
    ${LIB_PAIR_PREFIX}ptzblas
)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(ptzblas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

function(add_ptzblas_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # PTZBLAS depends on BLAS only — no BLACS/MPI cycle. We still wrap
    # in LINK_GROUP(RESCAN) because ptzblas/blas/refblas_quad form a
    # mutually-dependent static-archive set (refblas_quad provides the
    # quad-promoted DASUM/DDOT, ptzblas provides QVASUM, etc.).
    set(_ptzblas_cycle
        ptzblas_test_target ptzblas_test_common
        ${LIB_PAIR_PREFIX}ptzblas ${LIB_PAIR_PREFIX}blas refblas_quad)
    foreach(_common blas_common ptzblas_common)
        if(TARGET ${_common})
            list(APPEND _ptzblas_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_ptzblas_cycle}>")
    add_test(NAME ptzblas_${test_name}
             COMMAND $<TARGET_FILE:${test_name}>
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

file(GLOB _tests CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.f90)
foreach(_t ${_tests})
    add_ptzblas_precision_test(${_t})
endforeach()
