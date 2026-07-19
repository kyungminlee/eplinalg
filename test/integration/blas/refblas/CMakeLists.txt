# refblas_quad — quad-precision Netlib BLAS reference.
#
# Built from the vendored Netlib BLAS sources (copied to
# ${staging}/_refblas_src by `migrator stage`) compiled with gfortran's
# -freal-8-real-16, which promotes REAL(KIND=8) and DOUBLE PRECISION
# entities to REAL(KIND=16) without renaming routines. The resulting
# library exposes the standard D/Z routine names (DGEMM, DDOT, ZAXPY,
# ...) operating at quad precision — used as the universal reference
# for the differential precision tests.
#
# Selection (from parent tests/blas/CMakeLists.txt):
#   -DBLAS_REF=netlib  (default) → build from vendored sources
#   -DBLAS_REF=system            → link system -lblas (double, capped)

if(BLAS_REF STREQUAL "system")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

# Vendored-source mode requires gfortran for the -freal-8-real-16 flag.
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "tests/blas/refblas: BLAS_REF=netlib needs gfortran "
                   "for -freal-8-real-16 promotion. Falling back to "
                   "system -lblas.")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

set(_refsrc_dir ${CMAKE_SOURCE_DIR}/_refblas_src)
if(NOT IS_DIRECTORY ${_refsrc_dir})
    message(STATUS "tests/blas/refblas: ${_refsrc_dir} not found "
                   "(expected to be copied by `migrator stage`); "
                   "falling back to system -lblas.")
    add_library(refblas_quad INTERFACE)
    target_link_libraries(refblas_quad INTERFACE blas)
    return()
endif()

file(GLOB _refsrcs_f   CONFIGURE_DEPENDS ${_refsrc_dir}/*.f)
file(GLOB _refsrcs_f90 CONFIGURE_DEPENDS ${_refsrc_dir}/*.f90)
set(_refsrcs ${_refsrcs_f} ${_refsrcs_f90})
# Skip mixed-precision routines (single+double): not promotable cleanly.
list(FILTER _refsrcs EXCLUDE REGEX "/(dsdot|sdsdot)\\.f9?0?$")

add_library(refblas_quad STATIC ${_refsrcs})
# -freal-8-real-16 promotes REAL(KIND=8) and DOUBLE PRECISION (= KIND=8
# by default) entities to REAL(KIND=16), as well as their COMPLEX
# counterparts. D/Z literals (1.0D0, (1.0D0,2.0D0)) get promoted too.
# S/C single-precision routines are unaffected (KIND=4 stays KIND=4).
# -std=legacy is gfortran-only; -w silences the avalanche of "obsolete"
# warnings from the upstream f77 sources.
target_compile_options(refblas_quad PRIVATE
    -freal-8-real-16
    -w
    -std=legacy
)

# Rename Fortran-mangled exports from `<name>_` to `<name>_quad_` so
# refblas_quad doesn't collide with the standard ``blas`` archive's
# native-precision symbols when the kind4 / kind8 baseline tests link
# both. The matching update on the consumer side lives in
# tests/common/ref_quad_blas.f90: every routine declared in the
# interface block is suffixed ``_quad``, with a thin module procedure
# under the original name forwarding to the renamed external. Test
# files keep ``use ref_quad_blas, only: dasum`` unchanged.
#
# Looks for the rename script in two places — the in-tree path
# scripts/refquad_rename_archive.sh (when this CMakeLists is consumed
# from the repo) and ../scripts/... (when consumed from the staging
# tree, where ``migrator stage`` plants the cmake/ files alongside
# ``scripts/``). Skipped silently if neither path resolves.
find_file(REFQUAD_RENAME_SCRIPT refquad_rename_archive.sh
    PATHS
        ${CMAKE_SOURCE_DIR}/scripts
        ${CMAKE_SOURCE_DIR}/../scripts
        ${PROJECT_SOURCE_DIR}/scripts
        ${PROJECT_SOURCE_DIR}/../scripts
    NO_DEFAULT_PATH)
if(REFQUAD_RENAME_SCRIPT)
    # Common-archive excludes: their T/W defs are migrated-side
    # helpers (la_constants_ep, iparam2stage_ep_, mumps_abort_,
    # …). reflapack_quad references some of them; we must NOT rename
    # those references because the migrated link cycle expects them
    # under their original names.
    set(_refquad_excludes)
    foreach(_c blas_common lapack_common scalapack_common
               pblas_common pbblas_common ptzblas_common blacs_common
               mumps_common xblas_common)
        if(TARGET ${_c})
            list(APPEND _refquad_excludes -x $<TARGET_FILE:${_c}>)
        endif()
    endforeach()
    add_custom_command(TARGET refblas_quad POST_BUILD
        COMMAND sh ${REFQUAD_RENAME_SCRIPT}
                $<TARGET_FILE:refblas_quad>
                ${_refquad_excludes}
        COMMENT "refblas_quad: rename Fortran exports -> *_quad_"
        VERBATIM)
else()
    message(WARNING "refblas_quad: rename script not found; "
                    "kind4/kind8 baseline tests will see symbol collisions")
endif()
# .f90 files are free-form by extension; the rest are fixed-form.
# No per-file form override needed.
