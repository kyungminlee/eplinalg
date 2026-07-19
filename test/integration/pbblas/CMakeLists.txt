# Differential precision tests for migrated PBBLAS.
#
# PBBLAS is the distributed PB-BLAS auxiliary helper layer used inside
# PBLAS — Fortran with BLACS for the few routines that move data
# between processes (pbdtran/pbdtrnv) and pure local arithmetic for
# the rest (pbdmatadd/pbdvecadd).
#
# Each test program builds a small problem at REAL(KIND=16), runs on a
# BLACS grid (default 2x2 → mpiexec -n 4), calls the migrated
# q/e/dd-prefix entry point through a per-target wrapper, and either
#   - verifies locally on each rank (matadd/vecadd: data is local), or
#   - gathers the distributed result to rank 0 and compares against a
#     serial reference (pbdtran).

if(NOT HAVE_REAL16)
    message(STATUS "tests/pbblas: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/pbblas: skipping — MPI Fortran not found "
                   "(required for BLACS grid init)")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/pbblas: skipping — refblas_quad target not built "
                   "(tests/blas must be staged alongside tests/pbblas)")
    return()
endif()

if(NOT TARGET ${LIB_PAIR_PREFIX}pbblas)
    message(STATUS "tests/pbblas: skipping — ${LIB_PAIR_PREFIX}pbblas not built")
    return()
endif()

if(NOT TARGET ${LIB_PAIR_PREFIX}blacs)
    message(STATUS "tests/pbblas: skipping — ${LIB_PAIR_PREFIX}blacs not built")
    return()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# The shared tests/common helper modules are compiled once per suite
# (with per-suite reporter macros), so a per-tree mod dir keeps each
# tree's .mod files isolated in the shared fmod_tests/ tree.
set(_pbb_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}/pbblas_tests")
file(MAKE_DIRECTORY ${_pbb_moddir})

# ── Common helper modules ───────────────────────────────────────────
# pblas_distrib is a fypp template (one body instantiated for the real
# and complex families) — expand it into this suite's binary dir.
add_fypp_source(_pbb_distrib_src
    ${CMAKE_CURRENT_SOURCE_DIR}/../common/pblas_distrib.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/../common)
add_library(pbblas_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    ../common/prec_report.F90
    ../common/pbblas_prec_report.f90
    ../common/test_data.f90
    ../common/ref_quad_blas.f90
    ../common/pblas_grid.F90
    ${_pbb_distrib_src}
)
set_target_properties(pbblas_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_pbb_moddir}"
)
# MPI reporter: rank 0 writes the JSON, report_finalize broadcasts the
# pass/fail flag — see tests/common/prec_report.F90.
target_compile_definitions(pbblas_test_common PRIVATE PREC_REPORT_MPI)
target_include_directories(pbblas_test_common PUBLIC
    $<BUILD_INTERFACE:${_pbb_moddir}>
)
target_link_libraries(pbblas_test_common PUBLIC
    refblas_quad
    ${LIB_PAIR_PREFIX}blacs
    MPI::MPI_Fortran
)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/pbblas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_pbblas_wrapper_src
    ${_wrapper_dir}/target_pbblas.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(pbblas_test_target STATIC ${_pbblas_wrapper_src})
set_target_properties(pbblas_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_pbb_moddir}"
)
target_include_directories(pbblas_test_target PUBLIC
    $<BUILD_INTERFACE:${_pbb_moddir}>
)
target_link_libraries(pbblas_test_target PUBLIC
    pbblas_test_common
    target_conv
    ${LIB_PAIR_PREFIX}pbblas
)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(pbblas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# Default 4 ranks → 2×2 grid. Override at configure time:
# -DPBBLAS_TEST_NPROC=<count>.
set(PBBLAS_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/pbblas (default 4 → 2x2 grid)")

# Tests whose filename ends in `_lcm6` need a 6-rank world to form a
# 2×3 / 3×2 grid where MIN(LCMQ,CEIL(M,NB)) > 1 in pb[dz]tran.f. They
# self-skip via grid_init_shape's nproc-mismatch path when launched
# with a different rank count, so this knob is only load-bearing on
# presets that can deliver a connected MPI world (e.g. linux-impi).
set(PBBLAS_TEST_NPROC_LCM6 6 CACHE STRING
    "MPI process count for tests/pbblas *_lcm6 tests (default 6 → 2x3/3x2 grid)")

function(add_pbblas_precision_test test_src nproc)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # Same -fno-lto rationale as tests/pblas: Ubuntu's default
    # -flto=auto confuses the linker about which archive members
    # satisfy which cross-archive references.
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # PBBLAS / PTZBLAS / BLACS / BLAS form a cyclic static-archive
    # graph — pbblas references blacs (Cblacs_gridinfo / DGEBR2D-type
    # comm primitives) and ptzblas/blas; ptzblas references blas.
    # Wrap in LINK_GROUP(RESCAN) so the linker re-scans member
    # archives until cross-references settle.
    set(_pbblas_cycle
        pbblas_test_target pbblas_test_common
        ${LIB_PAIR_PREFIX}pbblas ${LIB_PAIR_PREFIX}ptzblas
        ${LIB_PAIR_PREFIX}blacs ${LIB_PAIR_PREFIX}blas refblas_quad)
    foreach(_common blas_common blacs_common ptzblas_common pbblas_common)
        if(TARGET ${_common})
            list(APPEND _pbblas_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_pbblas_cycle}>")
    add_test(NAME pbblas_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

file(GLOB _tests CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.f90)
foreach(_t ${_tests})
    get_filename_component(_tn ${_t} NAME_WE)
    if(_tn MATCHES "_lcm6$")
        add_pbblas_precision_test(${_t} ${PBBLAS_TEST_NPROC_LCM6})
    else()
        add_pbblas_precision_test(${_t} ${PBBLAS_TEST_NPROC})
    endif()
endforeach()
