# Differential precision tests for migrated ScaLAPACK.
#
# Mirrors tests/pblas/CMakeLists.txt: each test program generates a
# global problem at REAL(KIND=16), distributes block-cyclicly over a
# BLACS grid, calls the migrated q/x-prefix ScaLAPACK routine (via the
# per-target wrapper), gathers the distributed result to rank 0, and
# compares against serial reflapack_quad (or refblas_quad for pure-BLAS
# auxiliaries) applied to the full global matrix.
#
# ScaLAPACK depends on BLAS + LAPACK + BLACS + PBLAS, so all four
# migrated libraries must be present alongside the reference quad
# libraries for this suite to compile.

if(NOT HAVE_REAL16)
    message(STATUS "tests/scalapack: skipping — HAVE_REAL16 is required")
    return()
endif()

if(NOT MPI_Fortran_FOUND)
    message(STATUS "tests/scalapack: skipping — MPI Fortran not found")
    return()
endif()

if(NOT TARGET refblas_quad)
    message(STATUS "tests/scalapack: skipping — refblas_quad not available "
                   "(tests/blas must be staged)")
    return()
endif()

if(NOT TARGET reflapack_quad)
    message(STATUS "tests/scalapack: skipping — reflapack_quad not available "
                   "(tests/lapack must be staged)")
    return()
endif()

if(BASELINE_BUILD)
    # kind4/kind8: standard archives are the libraries under test.
    foreach(_lib scalapack pblas blacs lapack)
        if(NOT TARGET ${_lib})
            message(STATUS "tests/scalapack: skipping — ${_lib} std archive not built")
            return()
        endif()
    endforeach()
else()
    foreach(_lib ${LIB_PAIR_PREFIX}scalapack ${LIB_PAIR_PREFIX}pblas
                 ${LIB_PAIR_PREFIX}blacs ${LIB_PAIR_PREFIX}lapack)
        if(NOT TARGET ${_lib})
            message(STATUS "tests/scalapack: skipping — ${_lib} not built")
            return()
        endif()
    endforeach()
endif()

# ── Dedicated module directory ──────────────────────────────────────
# tests/scalapack/common reuses the same module names as
# tests/pblas/common (prec_kinds, compare, pblas_prec_report,
# test_data, ref_quad_blas, pblas_grid, pblas_distrib). Placing the
# ScaLAPACK variants' .mod files in their own subdirectory under fmod_tests
# prevents collisions when both test trees build together.
set(_sp_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}/scalapack_tests")
file(MAKE_DIRECTORY ${_sp_moddir})

# ── Common helper modules ───────────────────────────────────────────
add_library(scalapack_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    ../common/test_data.f90
    common/ref_quad_blas.f90
    common/ref_quad_lapack.f90
    common/cond_helpers.f90
    common/pblas_grid.f90
    ../common/pblas_distrib.f90
)
set_target_properties(scalapack_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_sp_moddir}"
)
target_include_directories(scalapack_test_common PUBLIC
    $<BUILD_INTERFACE:${_sp_moddir}>
)
if(BASELINE_BUILD)
    target_link_libraries(scalapack_test_common PUBLIC
        refblas_quad
        reflapack_quad
        blacs
        MPI::MPI_Fortran
    )
else()
    target_link_libraries(scalapack_test_common PUBLIC
        refblas_quad
        reflapack_quad
        ${LIB_PAIR_PREFIX}blacs
        MPI::MPI_Fortran
    )
endif()

# ── Per-target wrapper module ───────────────────────────────────────
set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/scalapack: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

add_fypp_source(_scalapack_wrapper_src
    ${_wrapper_dir}/target_scalapack.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(scalapack_test_target STATIC ${_scalapack_wrapper_src})
set_target_properties(scalapack_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_sp_moddir}"
)
target_include_directories(scalapack_test_target PUBLIC
    $<BUILD_INTERFACE:${_sp_moddir}>
)
if(BASELINE_BUILD)
    target_link_libraries(scalapack_test_target PUBLIC
        scalapack_test_common
        target_conv
        scalapack
        pblas
        lapack
        blas
    )
    foreach(_dep ptzblas pbblas)
        if(TARGET ${_dep})
            target_link_libraries(scalapack_test_target PUBLIC ${_dep})
        endif()
    endforeach()
else()
    target_link_libraries(scalapack_test_target PUBLIC
        scalapack_test_common
        target_conv
        ${LIB_PAIR_PREFIX}scalapack
        ${LIB_PAIR_PREFIX}pblas
        ${LIB_PAIR_PREFIX}lapack
        ${LIB_PAIR_PREFIX}blas
    )
    # qxscalapack's Fortran routines call qlamov_/pqgemr2d_ etc.; those C
    # helpers are now folded into ${LIB_PAIR_PREFIX}scalapack itself, so no
    # separate qscalapack_c link is needed.
    foreach(_dep ${LIB_PAIR_PREFIX}ptzblas ${LIB_PAIR_PREFIX}pbblas)
        if(TARGET ${_dep})
            target_link_libraries(scalapack_test_target PUBLIC ${_dep})
        endif()
    endforeach()
endif()
if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
    target_link_libraries(scalapack_test_target PUBLIC multifloatsf)
endif()

# ── Test executables ────────────────────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

set(SCALAPACK_TEST_NPROC 4 CACHE STRING
    "MPI process count for tests/scalapack (default 4 → 2x2 grid)")

function(add_scalapack_precision_test test_src)
    get_filename_component(test_name ${test_src} NAME_WE)
    add_executable(${test_name} ${test_src})
    # -fno-lto: same LTO-confusion issue documented in tests/pblas.
    target_compile_options(${test_name} PRIVATE -fno-lto)
    target_link_options(${test_name} PRIVATE -fno-lto)
    # ScaLAPACK ↔ LAPACK ↔ PBLAS ↔ PTZBLAS ↔ BLACS ↔ BLAS form a
    # cyclic archive graph. Wrap the whole cycle in LINK_GROUP(RESCAN)
    # so ld re-scans member archives until all cross-references resolve.
    if(BASELINE_BUILD)
        set(_sp_cycle
            scalapack_test_target scalapack_test_common
            scalapack pblas ptzblas pbblas
            lapack blas blacs
            reflapack_quad refblas_quad)
    else()
        set(_sp_cycle
            scalapack_test_target scalapack_test_common
            ${LIB_PAIR_PREFIX}scalapack ${LIB_PAIR_PREFIX}pblas
            ${LIB_PAIR_PREFIX}ptzblas  ${LIB_PAIR_PREFIX}pbblas
            ${LIB_PAIR_PREFIX}lapack   ${LIB_PAIR_PREFIX}blas
            ${LIB_PAIR_PREFIX}blacs
            reflapack_quad refblas_quad)
        # The dlamov-family wrappers (qlamov_, xlamov_) and REDIST entry
        # points (pqgemr2d_, pxgemr2d_) that ${LIB_PAIR_PREFIX}scalapack's
        # Fortran routines call are now folded into ${LIB_PAIR_PREFIX}scalapack
        # itself — no separate qscalapack_c archive to add.
        # qxxblas is a PUBLIC dep of qxlapack (cmake/CMakeLists.txt:600) for
        # the iterative-refinement extra-precision routines. It pulls qxblas,
        # which is already in the cycle — must include qxxblas explicitly so
        # the LINK_GROUP closes; otherwise CMake reports a strongly-connected
        # component cycle (group ↔ qxxblas).
        if(TARGET ${LIB_PAIR_PREFIX}xblas)
            list(APPEND _sp_cycle ${LIB_PAIR_PREFIX}xblas)
        endif()
    endif()
    foreach(_common blas_common blacs_common ptzblas_common pblas_common lapack_common scalapack_common)
        if(TARGET ${_common})
            list(APPEND _sp_cycle ${_common})
        endif()
    endforeach()
    target_link_libraries(${test_name} PRIVATE
        "$<LINK_GROUP:RESCAN,${_sp_cycle}>")
    add_test(NAME scalapack_${test_name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${SCALAPACK_TEST_NPROC}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${test_name}> ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
endfunction()

foreach(_group linear_solve factorization eigenvalue svd auxiliary)
    file(GLOB _tests CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/${_group}/test_*.f90)
    foreach(_t ${_tests})
        add_scalapack_precision_test(${_t})
    endforeach()
endforeach()
