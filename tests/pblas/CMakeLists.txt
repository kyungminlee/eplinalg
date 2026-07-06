# Differential precision tests for migrated PBLAS.
#
# Each test program generates a global problem at REAL(KIND=16),
# distributes block-cyclicly over a BLACS grid, calls the migrated
# q/e/dd-prefix PBLAS routine (via the per-target wrapper), gathers
# the distributed result to rank 0, and compares against serial
# refblas_quad applied to the full global matrix.
#
# The serial BLAS reference is supplied by the refblas_quad target
# built under tests/blas/refblas — i.e. tests/blas must also be
# staged for the PBLAS tests to compile.

if(NOT HAVE_REAL16)
    message(STATUS "tests/pblas: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/pblas: skipping — MPI Fortran not found "
                   "(required for BLACS grid init and rank-0 gather)")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/pblas: skipping — refblas_quad target not built "
                   "(tests/blas must be staged alongside tests/pblas)")
    return()
endif()

if(BASELINE_BUILD)
    # kind4/kind8: link standard Netlib pblas/blacs archives. Skip the
    # migrated-target presence checks since LIB_PAIR_PREFIX is empty.
    if(NOT TARGET pblas)
        message(STATUS "tests/pblas: skipping — pblas std archive not built")
        return()
    endif()
    if(NOT TARGET blacs)
        message(STATUS "tests/pblas: skipping — blacs std archive not built")
        return()
    endif()
else()
    if(NOT TARGET ${LIB_PAIR_PREFIX}pblas)
        message(STATUS "tests/pblas: skipping — ${LIB_PAIR_PREFIX}pblas not built")
        return()
    endif()
    if(NOT TARGET ${LIB_PAIR_PREFIX}blacs)
        message(STATUS "tests/pblas: skipping — ${LIB_PAIR_PREFIX}blacs not built")
        return()
    endif()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# pblas_grid / pblas_distrib module names collide with the (now
# divergent) copies in tests/pbblas and tests/blacs. Per-tree mod dir
# keeps each tree's .mod files isolated in the shared fmod_tests/ tree.
set(_pb_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}/pblas_tests")
file(MAKE_DIRECTORY ${_pb_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(pblas_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    ../common/test_data.f90
    common/ref_quad_blas.f90
    common/pblas_grid.f90
    ../common/pblas_distrib.f90
)
set_target_properties(pblas_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_pb_moddir}"
)
target_include_directories(pblas_test_common PUBLIC
    $<BUILD_INTERFACE:${_pb_moddir}>
)
if(BASELINE_BUILD)
    target_link_libraries(pblas_test_common PUBLIC
        refblas_quad
        blacs
        MPI::MPI_Fortran
    )
else()
    target_link_libraries(pblas_test_common PUBLIC
        refblas_quad
        ${LIB_PAIR_PREFIX}blacs
        MPI::MPI_Fortran
    )
endif()

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/pblas: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_pblas_wrapper_src
    ${_wrapper_dir}/target_pblas.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(pblas_test_target STATIC ${_pblas_wrapper_src})
set_target_properties(pblas_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_pb_moddir}"
)
target_include_directories(pblas_test_target PUBLIC
    $<BUILD_INTERFACE:${_pb_moddir}>
)
if(BASELINE_BUILD)
    target_link_libraries(pblas_test_target PUBLIC
        pblas_test_common
        target_conv
        pblas
    )
else()
    target_link_libraries(pblas_test_target PUBLIC
        pblas_test_common
        target_conv
        ${LIB_PAIR_PREFIX}pblas
    )
endif()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(pblas_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# Default 4 ranks → 2×2 grid (even for NP=1 this degenerates to 1×1
# and every test still passes, just as serial BLAS). Override at
# configure time: -DPBLAS_TEST_NPROC=<count>.
set(PBLAS_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/pblas (default 4 → 2x2 grid)")

function(add_pblas_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # Turn off LTO for the final link. Ubuntu's default -flto=auto
    # confuses the linker about which archive members satisfy which
    # cross-archive references (we've seen Cblacs_pnum / dwalltime00_
    # reported undefined even when --start-group is in effect). Fat
    # LTO objects still contain native code, so disabling LTO only
    # changes which representation the linker chooses.
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # PBLAS / PTZBLAS / BLACS / C-API shim form a cyclic static-archive
    # graph — e.g. libqpblas references Cblacs_gridinfo (in the shim),
    # the shim references MPI, PTZBLAS references BLAS/BLACS, PBLAS
    # references PTZBLAS. GNU ld in its default single-pass mode
    # leaves some symbols unresolved; LINK_GROUP(RESCAN) wraps the
    # set in --start-group/--end-group so the linker re-scans every
    # member archive until all cross-references settle.
    # LINK_GROUP(RESCAN, ...) only wraps the explicitly-named targets in
    # --start-group/--end-group — transitive deps end up outside the
    # group. We therefore enumerate every archive that participates in
    # the BLACS↔PBLAS↔PTZBLAS dependency cycle.
    if(BASELINE_BUILD)
        # `scalapack` carries TOOLS routines (zzdotc_, zzdotu_, …) that
        # PTZBLAS forwards calls into; PBLAS tests need it on the link
        # line even though they don't otherwise touch ScaLAPACK.
        set(_pblas_cycle
            pblas_test_target pblas_test_common
            pblas ptzblas pbblas
            blacs blas scalapack refblas_quad)
    else()
        set(_pblas_cycle
            pblas_test_target pblas_test_common
            ${LIB_PAIR_PREFIX}pblas ${LIB_PAIR_PREFIX}ptzblas ${LIB_PAIR_PREFIX}pbblas
            ${LIB_PAIR_PREFIX}blacs ${LIB_PAIR_PREFIX}blas refblas_quad)
    endif()
    # pblas_common / blacs_common / ptzblas_common hold the common
    # (precision-independent) halves split out by _install_library_pair;
    # they are referenced inside the cycle too.
    foreach(_common blas_common blacs_common ptzblas_common pblas_common)
        if(TARGET ${_common})
            list(APPEND _pblas_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_pblas_cycle}>")
    add_test(NAME pblas_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PBLAS_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_lvl level1 level2 level3)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_lvl}/test_*.f90)
    foreach(_t ${_tests})
        add_pblas_precision_test(${_t})
    endforeach()
endforeach()
