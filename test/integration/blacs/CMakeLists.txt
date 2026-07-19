# Differential precision tests for migrated BLACS.
#
# BLACS is a thin Fortran-callable layer over MPI providing process
# grid management plus per-type point-to-point send/recv, broadcast,
# and combine (sum/amx/amn) primitives. Communication is exact; we
# verify byte-equivalent payloads at REAL(KIND=16) using a 0.0
# tolerance (no roundoff).
#
# Each test program runs on a 2x2 BLACS grid (default mpiexec -n 4),
# prepares known per-rank data, calls the migrated entry point through
# a per-target wrapper, and asserts exact equality.

if(NOT HAVE_REAL16)
    message(STATUS "tests/blacs: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/blacs: skipping — MPI Fortran not found")
    return()
endif()

if(NOT TARGET ${LIB_PAIR_PREFIX}blacs)
    message(STATUS "tests/blacs: skipping — ${LIB_PAIR_PREFIX}blacs not built")
    return()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# The shared tests/common helper modules are compiled once per suite
# (with per-suite reporter macros), so a per-tree subdir keeps our
# .mod files from colliding in the shared fmod_tests/ tree.
set(_bl_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}/blacs_tests")
file(MAKE_DIRECTORY ${_bl_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(blacs_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    ../common/prec_report.F90
    ../common/blacs_prec_report.f90
    ../common/test_data.f90
    ../common/pblas_grid.F90
)
set_target_properties(blacs_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_bl_moddir}"
)
# MPI reporter: rank 0 writes the JSON, report_finalize broadcasts the
# pass/fail flag — see tests/common/prec_report.F90.
target_compile_definitions(blacs_test_common PRIVATE PREC_REPORT_MPI)
target_include_directories(blacs_test_common PUBLIC
    $<BUILD_INTERFACE:${_bl_moddir}>
)
target_link_libraries(blacs_test_common PUBLIC
    ${LIB_PAIR_PREFIX}blacs
    MPI::MPI_Fortran
)

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/blacs: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_blacs_wrapper_src
    ${_wrapper_dir}/target_blacs.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(blacs_test_target STATIC ${_blacs_wrapper_src})
set_target_properties(blacs_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_bl_moddir}"
)
target_include_directories(blacs_test_target PUBLIC
    $<BUILD_INTERFACE:${_bl_moddir}>
)
target_link_libraries(blacs_test_target PUBLIC
    blacs_test_common
    target_conv
    ${LIB_PAIR_PREFIX}blacs
)
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(blacs_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# Default 4 ranks → 2x2 grid. Override at configure time:
# -DBLACS_TEST_NPROC=<count>.
set(BLACS_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/blacs (default 4 → 2x2 grid)")

function(add_blacs_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    set(_blacs_cycle
        blacs_test_target blacs_test_common
        ${LIB_PAIR_PREFIX}blacs)
    foreach(_common blas_common blacs_common)
        if(TARGET ${_common})
            list(APPEND _blacs_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_blacs_cycle}>")
    add_test(NAME blacs_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${BLACS_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

file(GLOB _tests CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.F90)
foreach(_t ${_tests})
    add_blacs_precision_test(${_t})
endforeach()
