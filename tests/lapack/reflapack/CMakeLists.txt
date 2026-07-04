# reflapack_quad — quad-precision Netlib LAPACK reference.
#
# Built from the vendored Netlib LAPACK sources (copied to
# ${staging}/_reflapack_src by `migrator stage`) compiled with gfortran's
# -freal-8-real-16, which promotes REAL(KIND=8) and DOUBLE PRECISION
# entities to REAL(KIND=16) without renaming routines. The resulting
# library exposes the standard D/Z routine names (DGESV, DGETRF, ZHEEV,
# ...) operating at quad precision — used as the universal reference
# for the differential precision tests.
#
# Selection (from parent tests/lapack/CMakeLists.txt):
#   -DLAPACK_REF=netlib  (default) → build from vendored sources
#   -DLAPACK_REF=system            → link system -llapack (double, capped)

if(LAPACK_REF STREQUAL "system")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

# Vendored-source mode requires gfortran for the -freal-8-real-16 flag.
if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    message(STATUS "tests/lapack/reflapack: LAPACK_REF=netlib needs gfortran "
                   "for -freal-8-real-16 promotion. Falling back to "
                   "system -llapack.")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

set(_reflap_dir ${CMAKE_SOURCE_DIR}/_reflapack_src)
if(NOT IS_DIRECTORY ${_reflap_dir})
    message(STATUS "tests/lapack/reflapack: ${_reflap_dir} not found "
                   "(expected to be copied by `migrator stage`); "
                   "falling back to system -llapack.")
    add_library(reflapack_quad INTERFACE)
    target_link_libraries(reflapack_quad INTERFACE lapack refblas_quad)
    return()
endif()

file(GLOB _reflap_f   CONFIGURE_DEPENDS ${_reflap_dir}/*.f)
file(GLOB _reflap_F   CONFIGURE_DEPENDS ${_reflap_dir}/*.F)
file(GLOB _reflap_f90 CONFIGURE_DEPENDS ${_reflap_dir}/*.f90)
file(GLOB _reflap_F90 CONFIGURE_DEPENDS ${_reflap_dir}/*.F90)
set(_reflap_srcs ${_reflap_f} ${_reflap_F} ${_reflap_f90} ${_reflap_F90})

# Skip mixed-precision / cross-precision routines (single+double or
# complex+real32 bridges) that don't promote cleanly under
# -freal-8-real-16. Match either case of the .f / .F suffix in case
# upstream Netlib ever renames a bridge file to capital F.
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dsgesv|zcgesv|dsposv|zcposv|dsgels|zcgels)\\.[fF]9?0?$")
list(FILTER _reflap_srcs EXCLUDE REGEX "/(dlag2s|slag2d|zlag2c|clag2z|dlat2s|zlat2c)\\.[fF]9?0?$")

# Skip multifloats-specific auxiliary modules (la_*_mf / la_*_ep are
# for targets that don't use the standard -freal-8-real-16 promotion).
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_constants_ep")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_mf")
list(FILTER _reflap_srcs EXCLUDE REGEX "la_xisnan_ep")

# Skip iparam2stage.F. It is a precision-independent INTEGER function
# (no FP argument promotion), and the staged lapack/src/iparam2stage.F
# already ships a patched copy via recipes/lapack/mf_overrides/ that
# recognizes the migrator precision prefixes (Q / X / E / Y / M / W).
# If we let reflapack_quad keep its own iparam2stage_ symbol the linker
# resolves it from this archive first, which masks the patched copy in
# liblapack_common — and qsyev_2stage / xheev_2stage / etc. then get
# IPARAM2STAGE = -1 on the migrated path, undersized LWORK, and heap
# corruption inside libc free. Excluding the upstream copy here lets
# the linker pull the patched version from liblapack_common for both
# the reference (D-prefix) and migrated (Q/E/M/X/Y/W-prefix) call
# sites — they're ABI-compatible because iparam2stage takes only
# INTEGER and CHARACTER arguments.
list(FILTER _reflap_srcs EXCLUDE REGEX "/iparam2stage\\.F$")

# Reference-source overrides. Files in overrides/ replace their
# same-named counterparts from the upstream Netlib SRC tree. Used to
# patch latent precision bugs in the reference path that surface only
# under -freal-8-real-16 promotion (see overrides/zgedmd.f90 for the
# CMPLX(...,KIND=WP) fp64 truncation that broke test_zgedmd /
# test_zgedmdq on kind16 — the migrator's xgedmd is correct; only the
# quad-promoted reference loses precision).
file(GLOB _ovr_srcs CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/overrides/*.f
    ${CMAKE_CURRENT_SOURCE_DIR}/overrides/*.F
    ${CMAKE_CURRENT_SOURCE_DIR}/overrides/*.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/overrides/*.F90)
foreach(_ovr ${_ovr_srcs})
    get_filename_component(_ovr_name ${_ovr} NAME)
    list(FILTER _reflap_srcs EXCLUDE REGEX "/${_ovr_name}$")
endforeach()
list(APPEND _reflap_srcs ${_ovr_srcs})

add_library(reflapack_quad STATIC ${_reflap_srcs})
target_compile_options(reflapack_quad PRIVATE
    -freal-8-real-16
    -w
    -std=legacy
)

# Rename exports `<name>_` -> `<name>_quad_` (mirrors refblas_quad —
# see tests/blas/refblas/CMakeLists.txt for rationale). Consumer side:
# tests/lapack/common/ref_quad_lapack.f90 (and the analogous PBLAS /
# ScaLAPACK ref_quad_*.f90) get a thin module-procedure forward layer.
find_file(REFQUAD_RENAME_SCRIPT refquad_rename_archive.sh
    PATHS
        ${CMAKE_SOURCE_DIR}/scripts
        ${CMAKE_SOURCE_DIR}/../scripts
        ${PROJECT_SOURCE_DIR}/scripts
        ${PROJECT_SOURCE_DIR}/../scripts
    NO_DEFAULT_PATH)
if(REFQUAD_RENAME_SCRIPT)
    # Pass refblas_quad as a sibling (its defs feed the rename set, so
    # reflapack_quad's U refs to dnrm2_ etc. get renamed in lockstep).
    # `-x` excludes flag the migrated `*_common` archives: their T/W
    # defs (la_constants_ep, iparam2stage_ patched, mumps_abort_, …)
    # are referenced by reflapack_quad and must stay un-renamed so the
    # migrated link cycle resolves them under their original names.
    set(_refquad_excludes)
    foreach(_c blas_common lapack_common scalapack_common
               pblas_common pbblas_common ptzblas_common blacs_common
               mumps_common xblas_common)
        if(TARGET ${_c})
            list(APPEND _refquad_excludes -x $<TARGET_FILE:${_c}>)
        endif()
    endforeach()
    add_custom_command(TARGET reflapack_quad POST_BUILD
        COMMAND sh ${REFQUAD_RENAME_SCRIPT}
                $<TARGET_FILE:reflapack_quad>
                $<TARGET_FILE:refblas_quad>
                ${_refquad_excludes}
        COMMENT "reflapack_quad: rename Fortran symbols -> *_quad_"
        VERBATIM)
    add_dependencies(reflapack_quad refblas_quad)
else()
    message(WARNING "reflapack_quad: rename script not found; "
                    "kind4/kind8 baseline tests will see symbol collisions")
endif()
# LAPACK internally calls BLAS. Link the quad-promoted refblas_quad so
# every routine in the reference chain operates at KIND=16.
target_link_libraries(reflapack_quad PUBLIC refblas_quad)

# xx-family iterative refinement (dgesvxx / zgesvxx / dgerfsx / ...) is
# implemented in ${_reflap_dir}/{d,z}la_*_extended.f. Those files call
# XBLAS Level-2 entry points by their upstream d/z names
# (blas_dgemv_x_, blas_dgbmv_x_, blas_dsymv_x_, blas_zhemv_x_, ...).
# We don't ship a standard-precision XBLAS archive, and the migrated
# ${LIB_PREFIX}xblas is at the *target* working precision (e.g. 80-bit
# extended on kind10, double-double on multifloats) — so its q/x/e/m/y
# entry points cannot satisfy reflapack_quad's d/z calls under
# -freal-8-real-16 (which promotes both to KIND=16 quad).
#
# Bridge with quad-precision Fortran shims in
# refxblas_quad_bridge.f90. The shims expose the standard
# blas_dgemv_x_ / blas_zhemv2_x_ / etc. names and forward to the
# quad-promoted refblas_quad routines (DGEMV/ZGEMV/DGBMV/ZGBMV/DSYMV/
# ZHEMV) for everything that has a BLAS analogue, hand-coding the
# complex-symmetric matvec which has none. At quad working precision
# XBLAS's "extra-precise" residual refinement is the identity, so the
# shim form (no head/tail decomposition) is mathematically faithful.
#
# Compiling the bridge with the same -freal-8-real-16 / KIND=16 ABI as
# reflapack_quad keeps the parameter-passing convention identical to
# the upstream LAPACK call sites.
target_sources(reflapack_quad PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/refxblas_quad_bridge.f90)

# kind4 / kind8 baseline: std lapack's iterative-refinement helpers
# (dla_*_extended.f, zla_*_extended.f) reference XBLAS entry points
# (blas_dgemv_x_, blas_dgbmv_x_, …). Those .o files get pulled into
# the test executable transitively via target_lapack's public
# wrappers (target_dgerfsx, target_dgesvxx, …) — the wrappers' xx /
# rfsx tests are filtered out of the baseline ctest set, but the
# wrappers themselves still live in the same target_lapack.f90.o.
# Provide aborting stubs in a separate archive so link succeeds; the
# stubs are unreachable at runtime under the baseline test mix.
if(BASELINE_BUILD)
    add_library(xblas_baseline_stubs STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/xblas_baseline_stubs.c)
    # --whole-archive forces ld to extract every stub .o regardless of
    # current undefined-symbol demand. Without it, ld scans the stub
    # archive before std lapack pulls in dla_*_extended.f.o (which
    # references blas_*_x_), so the stubs never get extracted and the
    # references stay unresolved. Plain PUBLIC propagation then falls
    # short — once ld passes the stubs in the link line it doesn't
    # re-scan them when std lapack later demands their symbols.
    target_link_libraries(reflapack_quad PUBLIC
        -Wl,--whole-archive xblas_baseline_stubs -Wl,--no-whole-archive)
endif()
