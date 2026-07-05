# Unified CMake build for migrated extended-precision libraries.
#
# This file is copied into the staging directory by:
#   python -m migrator stage <staging_dir> --target <target>
#
# Usage:
#   cmake -S <staging_dir> -B <staging_dir>/build -DCMAKE_BUILD_TYPE=Release
#   cmake --build <staging_dir>/build -j
#   cmake --install <staging_dir>/build --prefix /usr/local
#
cmake_minimum_required(VERSION 3.20)

# Load target configuration (written by the Python migration tool)
include(${CMAKE_CURRENT_SOURCE_DIR}/target_config.cmake)

# BASELINE_BUILD is true for the un-migrated kind4 / kind8 targets:
# `migrator stage --target kind4|kind8` writes target_config.cmake
# with LIB_PREFIX="" so the test framework links the standard archives
# (blas / lapack / pblas / scalapack) directly. We skip every
# add_migrated_* call and its post-link wiring under that mode — the
# bare archive names would otherwise collide with the standard archive
# of the same name (e.g. ``target_link_libraries(blas PUBLIC blas)``).
if("${LIB_PREFIX}" STREQUAL "")
    set(BASELINE_BUILD TRUE)
else()
    set(BASELINE_BUILD FALSE)
endif()

project(${TARGET_NAME}_libraries Fortran C CXX)

# ── Compiler-aware install infrastructure ────────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include(FortranCompiler)

# ── Global compiler flags ────────────────────────────────────────────
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
# Pin C standard to gnu17. The migrated ScaLAPACK C sources retain
# K&R-style empty parameter lists (``void f();``) which the C migrator
# treats as legal under ANSI C semantics. C23 (the default on GCC 15
# and recent Clang/icx) reinterprets ``()`` as ``(void)`` and turns
# call sites with arguments into hard errors, not warnings — so -w
# alone does not suppress them. Forcing gnu17 keeps the legacy
# semantics across GCC, Clang, icx, and flang.
set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS ON)
# The migrated C sources also call BLAS routines (e.g. zscal_) without
# explicit prototypes, relying on K&R implicit-int conventions. GCC 14+
# and recent Clang promote ``-Wimplicit-function-declaration`` to a
# hard error by default; downgrade it to a warning so the global -w
# flag can suppress it. Same story for ``-Wint-conversion`` (also an
# error by default in GCC 14+).
include(CheckCCompilerFlag)
foreach(_flag -Wno-error=implicit-function-declaration
              -Wno-error=int-conversion
              -Wno-error=incompatible-pointer-types)
    string(MAKE_C_IDENTIFIER "HAS_C${_flag}" _var)
    check_c_compiler_flag(${_flag} ${_var})
    if(${_var})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${_flag}")
    endif()
endforeach()
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
endif()
set(CMAKE_Fortran_PREPROCESS ON)

include(${CMAKE_CURRENT_SOURCE_DIR}/DetectExtendedPrecision.cmake)

find_package(MPI COMPONENTS C CXX Fortran QUIET)

# ── MPI flavor + version detection ───────────────────────────────────
# MPI-dependent libraries are not ABI-compatible across implementations
# (Intel MPI / OpenMPI / MPICH expose different ABIs for the same
# handles — see the libmpiseq comment below for the canonical example).
# Tag every MPI-dependent install with the implementation flavor and
# its major.minor version so multiple builds can coexist in the same
# install prefix and ``find_package(<lib>)`` can pick the matching
# variant.
#
# Probe mpi.h's vendor-specific version macros at configure time:
#   Intel MPI : I_MPI_VERSION (string, e.g. "2021.18.0")
#   OpenMPI   : OMPI_MAJOR_VERSION / OMPI_MINOR_VERSION (ints)
#   MPICH     : MPICH_VERSION (string, e.g. "4.2.0")
# Produces ``MPI_TAG`` like ``intelmpi-2021.18`` / ``openmpi-5.0`` /
# ``mpich-4.2``. Empty when no MPI was found.
set(MPI_TAG "")
if(MPI_C_FOUND)
    set(_mpi_inc ${MPI_C_HEADER_DIR})
    if(TARGET MPI::MPI_C)
        get_target_property(_iface MPI::MPI_C INTERFACE_INCLUDE_DIRECTORIES)
        if(_iface)
            list(APPEND _mpi_inc ${_iface})
        endif()
    endif()
    list(REMOVE_DUPLICATES _mpi_inc)
    set(_mpi_probe_c "${CMAKE_BINARY_DIR}/CMakeFiles/mpi_probe.c")
    file(WRITE "${_mpi_probe_c}" "\
#include <mpi.h>
#include <stdio.h>
int main(void) {
#if defined(I_MPI_VERSION)
    printf(\"intelmpi %s\\n\", I_MPI_VERSION);
#elif defined(OMPI_MAJOR_VERSION) && defined(OMPI_MINOR_VERSION)
    printf(\"openmpi %d.%d\\n\", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION);
#elif defined(MPICH_VERSION)
    printf(\"mpich %s\\n\", MPICH_VERSION);
#else
    printf(\"unknown ?\\n\");
#endif
    return 0;
}
")
    try_run(_mpi_run_rc _mpi_compile_rc
        "${CMAKE_BINARY_DIR}/CMakeFiles/mpi_probe.dir"
        "${_mpi_probe_c}"
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${_mpi_inc}"
        RUN_OUTPUT_VARIABLE _mpi_id)
    if(_mpi_compile_rc AND "${_mpi_run_rc}" STREQUAL "0")
        string(STRIP "${_mpi_id}" _mpi_id)
        # ``flavor major.minor[.patch]`` → ``flavor-major.minor``
        if(_mpi_id MATCHES "^([a-z]+) ([0-9]+)\\.([0-9]+)")
            set(MPI_TAG "${CMAKE_MATCH_1}-${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
        endif()
    endif()
    unset(_mpi_id)
    unset(_mpi_co)
endif()
if(MPI_TAG)
    message(STATUS "MPI tag: ${MPI_TAG}")
elseif(MPI_C_FOUND)
    message(STATUS "MPI found but could not identify flavor/version — MPI-dependent libraries will not be tag-suffixed")
endif()

# ── libmpiseq compatibility gate ─────────────────────────────────────
# libmpiseq's stubs (cmake/mpiseq_c_stubs.c et al.) are written against
# Intel MPI's handle ABI: MPI_BYTE, MPI_FLOAT, MPI_OP_NULL, etc. are
# integer constants in Intel MPI's mpi.h, so the stubs resolve cleanly
# at link time. OpenMPI and MPICH expose those handles as external
# symbol references (``ompi_mpi_byte``, ``mpich_mpi_byte``); linking
# libmpiseq.a against an executable on those MPIs produces unresolved
# references like ``undefined reference to `ompi_mpi_byte```.
#
# To prevent that silent failure mode (which previously surfaced as a
# build failure at the first ``_seq`` test executable), detect at
# configure time whether the active MPI is Intel MPI and only build
# the libmpiseq target on that host. Downstream test variants
# (``mumps_test_*_seq`` etc.) gate on ``TARGET mpiseq``, so they
# silently drop out on OpenMPI/MPICH hosts instead of failing to link.
if(MPI_TAG MATCHES "^intelmpi")
    set(MPISEQ_IS_INTEL_MPI TRUE)
else()
    set(MPISEQ_IS_INTEL_MPI FALSE)
endif()
if(MPISEQ_IS_INTEL_MPI)
    message(STATUS "libmpiseq: Intel MPI detected — _seq test variants enabled")
else()
    message(STATUS
        "libmpiseq: skipping (active MPI is not Intel MPI). "
        "_seq test variants will be omitted. Use the linux-impi CMake "
        "preset to enable them.")
endif()

# ── Multifloats (conditional) ────────────────────────────────────────
if(NEEDS_MULTIFLOATS)
    find_program(FYPP_EXECUTABLE fypp REQUIRED)

    # Use the multifloats top-level directory (not src/ and fsrc/
    # separately): upstream's src/CMakeLists.txt uses CMAKE_SOURCE_DIR
    # which is only correct when multifloats is the top-level project.
    # EXCLUDE_FROM_ALL keeps us from building multifloats's install
    # rules / tests by default.
    #
    # BUILD_TESTING/MULTIFLOATS_BUILD_BENCH are set as plain (non-cache)
    # variables so they only suppress tests *inside* the multifloats
    # sub-build. A `CACHE … FORCE` here would also disable our own
    # `tests/` directory at the end of this file.
    set(BUILD_TESTING OFF)
    set(MULTIFLOATS_BUILD_BENCH OFF)
    if(DEFINED MULTIFLOATS_DIR)
        message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
        add_subdirectory(${MULTIFLOATS_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
    else()
        include(FetchContent)
        set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
            CACHE STRING "Git URL for the multifloats library")
        # multifloats v0.6.0 fixed the ${CMAKE_SOURCE_DIR} include-path
        # leak that broke add_subdirectory'd builds (upstream issue #23).
        # Earlier tags need the file(READ ... REPLACE ... WRITE) shim
        # that used to live here; bumping below this floor is not
        # supported without restoring it.
        set(MULTIFLOATS_GIT_TAG "v0.6.0"
            CACHE STRING "Git tag/branch/commit for multifloats (>= v0.6.0)")
        message(STATUS "Fetching multifloats from ${MULTIFLOATS_GIT_REPO} (${MULTIFLOATS_GIT_TAG})")
        FetchContent_Declare(multifloats_fetch
            GIT_REPOSITORY ${MULTIFLOATS_GIT_REPO}
            GIT_TAG        ${MULTIFLOATS_GIT_TAG}
        )
        FetchContent_Populate(multifloats_fetch)
        add_subdirectory(
            ${multifloats_fetch_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
    endif()

    # la_constants_mw helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_constants_mw.f90")
        add_library(la_constants_mw STATIC "${MF_HELPERS_DIR}/la_constants_mw.f90")
        fortran_module_layout(la_constants_mw)
        target_link_libraries(la_constants_mw PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()

    # la_xisnan_mw helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_xisnan_mw.f90")
        add_library(la_xisnan_mw STATIC "${MF_HELPERS_DIR}/la_xisnan_mw.f90")
        fortran_module_layout(la_xisnan_mw)
        target_link_libraries(la_xisnan_mw PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()

    # multifloats_mpi — C++ MPI handle registration + bridge header
    if(EXISTS "${MF_HELPERS_DIR}/multifloats_mpi.cpp" AND MPI_CXX_FOUND)
        add_library(multifloats_mpi STATIC "${MF_HELPERS_DIR}/multifloats_mpi.cpp")
        target_include_directories(multifloats_mpi PUBLIC
            $<BUILD_INTERFACE:${MF_HELPERS_DIR}>)
        # MPI::MPI_CXX is PRIVATE so consumers (mblacs, mlapack PUBLIC
        # chain, mumps test exes) aren't forced to link libmpicxx.so.
        # Each consumer that wants the real-MPI runtime adds MPI::* on
        # its own — keeps the libmpiseq link path (which substitutes
        # mpiseq for ALL of MPI/BLACS/PBLAS/scalapack) clean.
        target_link_libraries(multifloats_mpi
            PRIVATE MPI::MPI_CXX
            PUBLIC  multifloats)
        target_compile_features(multifloats_mpi PRIVATE cxx_std_17)
    elseif(EXISTS "${MF_HELPERS_DIR}/multifloats_bridge.h")
        # No MPI — still provide the bridge header via an interface target
        add_library(multifloats_mpi INTERFACE)
        target_include_directories(multifloats_mpi INTERFACE
            $<BUILD_INTERFACE:${MF_HELPERS_DIR}>)
    endif()

    # multifloats_mpi_f — Fortran module exposing the C-side MPI_FLOAT64X2
    # / MPI_DD_SUM / etc. handles to Fortran callers (MUMPS calls MPI
    # from Fortran directly with these tokens; without the module they
    # have no IMPLICIT type at compile time). Linking into multifloats_mpi
    # ensures the C++ runtime that defines mf_mpi_*_f globals is pulled
    # in transitively whenever a target depends on multifloats_mpi_f.
    if(EXISTS "${MF_HELPERS_DIR}/multifloats_mpi_f.f90" AND TARGET multifloats_mpi)
        add_library(multifloats_mpi_f STATIC
            "${MF_HELPERS_DIR}/multifloats_mpi_f.f90")
        fortran_module_layout(multifloats_mpi_f)
        target_link_libraries(multifloats_mpi_f PUBLIC multifloats_mpi)
    endif()
endif()

# ── Quad MPI (conditional; kind16 REAL(KIND=16) / __float128) ─────────
# Custom MPI reduce ops (MPI_QQ_SUM / QQ_AMX / QQ_AMN and the XX_*
# complex analogues) registered on the *standard* MPI_REAL16 /
# MPI_COMPLEX32 datatypes. Intel MPI's builtin MPI_SUM/MAX/MIN have no
# 16-byte-real kernel and segfault on MPI_REDUCE at np>=2, so the
# migrated kind16 MUMPS routes every reduction through these instead —
# the quad analogue of the multifloats bridge, but plain C (no C++, no
# FetchContent, no custom datatypes).
if(NEEDS_QUAD_MPI)
    # quad_mpi — C registration of the custom combine ops. fabsq lives
    # in libquadmath.
    if(EXISTS "${MF_HELPERS_DIR}/quad_mpi.c")
        add_library(quad_mpi STATIC "${MF_HELPERS_DIR}/quad_mpi.c")
        target_include_directories(quad_mpi PUBLIC
            $<BUILD_INTERFACE:${MF_HELPERS_DIR}>)
        # MPI::MPI_C is PRIVATE so consumers linking libmpiseq instead of
        # the real MPI runtime (the _seq test variants) aren't forced to
        # pull in the real MPI .so — mirrors multifloats_mpi's handling.
        # (MPI_Op_c2f is a cast macro and MPI_Op_create is stubbed by
        # libmpiseq, so quad_mpi links cleanly against mpiseq too.)
        if(MPI_C_FOUND)
            target_link_libraries(quad_mpi PRIVATE MPI::MPI_C)
        endif()
        target_link_libraries(quad_mpi PRIVATE quadmath)
    endif()

    # quad_mpi_f — Fortran shim exposing the C-side op handles to MUMPS
    # (which issues MPI_REDUCE from Fortran with these tokens directly).
    # Linking quad_mpi PUBLIC pulls the C runtime that defines the *_f
    # globals into any target that depends on quad_mpi_f.
    if(EXISTS "${MF_HELPERS_DIR}/quad_mpi_f.f90" AND TARGET quad_mpi)
        add_library(quad_mpi_f STATIC "${MF_HELPERS_DIR}/quad_mpi_f.f90")
        fortran_module_layout(quad_mpi_f)
        target_link_libraries(quad_mpi_f PUBLIC quad_mpi)
    endif()
endif()

# ── Helper functions ─────────────────────────────────────────────────

# Build a stock, standard-precision archive (S/D/C/Z entry points) from
# vendored upstream sources. The target name is the bare library name
# (``blas``, ``lapack``, ...) — no ``_std`` suffix — so consumers can
# write ``target_link_libraries(... blas)`` and get the upstream
# DGEMM_/ZGEMM_/etc. symbols. The migrated archive (``${LIB_PREFIX}<name>``)
# PUBLIC-links it so Q/X/E/Y/M/W clones can call into the stock half
# without bundling the originals into the migrated archive.
#
# Args:
#   name      - Target name (bare, no prefix).
#   src_dir   - Directory containing upstream Fortran sources.
#   exclude_re - Optional regex passed to list(FILTER EXCLUDE) to drop
#                files that the migrator's recipe also drops (e.g.
#                cross-precision routines like dsdot/sdsdot).
function(add_standard_fortran_library name src_dir)
    cmake_parse_arguments(_std "" "EXCLUDE_REGEX" "" ${ARGN})
    if(NOT IS_DIRECTORY "${src_dir}")
        message(STATUS "${name} (standard): ${src_dir} not found — skipping")
        return()
    endif()
    # Glob both lower- and upper-case Fortran source extensions: ``.F``
    # and ``.F90`` are preprocessed Fortran (LAPACK uses capital-F for
    # iparam2stage.F and the la_*_ep / la_*_mf helpers). On case-
    # sensitive filesystems the lower-case glob misses them entirely.
    file(GLOB _f   CONFIGURE_DEPENDS ${src_dir}/*.f)
    file(GLOB _F   CONFIGURE_DEPENDS ${src_dir}/*.F)
    file(GLOB _f90 CONFIGURE_DEPENDS ${src_dir}/*.f90)
    file(GLOB _F90 CONFIGURE_DEPENDS ${src_dir}/*.F90)
    set(_srcs ${_f} ${_F} ${_f90} ${_F90})
    if(_std_EXCLUDE_REGEX)
        list(FILTER _srcs EXCLUDE REGEX "${_std_EXCLUDE_REGEX}")
    endif()
    if(NOT _srcs)
        message(STATUS "${name} (standard): no sources under ${src_dir} — skipping")
        return()
    endif()
    add_library(${name} STATIC ${_srcs})
    fortran_module_layout(${name})
endfunction()


# Build a stock, standard-precision C library archive (S/D/C/Z entry
# points) from vendored upstream sources. Sibling of
# ``add_standard_fortran_library`` for C libraries (BLACS, PBLAS).
# (ScaLAPACK's C-side sources are not built here — they are compiled
# as OBJECT libraries and folded directly into the ScaLAPACK archive;
# see the scalapack_c fold below.) Includes MPI privately and
# propagates it via INTERFACE for downstream Fortran link, mirroring
# add_migrated_c_library.
#
# Args:
#   name             - Target name (bare, no prefix).
#   src_dir          - Directory containing upstream C sources.
#   include_dirs     - Optional INCLUDE_DIRS list (defaults to src_dir).
#   compile_defs     - Optional COMPILE_DEFINITIONS list (e.g. ``Add_``).
#   exclude_re       - Optional EXCLUDE_REGEX passed to list(FILTER).
function(add_standard_c_library name src_dir)
    cmake_parse_arguments(_std "" "EXCLUDE_REGEX"
        "COMPILE_DEFINITIONS;INCLUDE_DIRS" ${ARGN})
    if(NOT IS_DIRECTORY "${src_dir}")
        message(STATUS "${name} (standard): ${src_dir} not found — skipping")
        return()
    endif()
    file(GLOB _srcs CONFIGURE_DEPENDS ${src_dir}/*.c)
    if(_std_EXCLUDE_REGEX)
        list(FILTER _srcs EXCLUDE REGEX "${_std_EXCLUDE_REGEX}")
    endif()
    if(NOT _srcs)
        message(STATUS "${name} (standard): no C sources under ${src_dir} — skipping")
        return()
    endif()
    add_library(${name} STATIC ${_srcs})
    if(_std_INCLUDE_DIRS)
        target_include_directories(${name} PUBLIC ${_std_INCLUDE_DIRS})
    else()
        target_include_directories(${name} PUBLIC $<BUILD_INTERFACE:${src_dir}>)
    endif()
    if(_std_COMPILE_DEFINITIONS)
        target_compile_definitions(${name} PRIVATE ${_std_COMPILE_DEFINITIONS})
    endif()
    if(MPI_C_FOUND)
        target_link_libraries(${name} PRIVATE MPI::MPI_C)
        target_link_libraries(${name} INTERFACE $<LINK_ONLY:MPI::MPI_C>)
    endif()
endfunction()


function(add_migrated_fortran_library lib_name)
    set(_dir ${CMAKE_CURRENT_SOURCE_DIR}/${lib_name})
    if(NOT EXISTS "${_dir}/manifest.cmake")
        return()
    endif()
    include(${_dir}/manifest.cmake)

    # Prefix source paths with library directory
    set(_common "")
    foreach(f ${${lib_name}_COMMON_SOURCES})
        list(APPEND _common "${_dir}/${f}")
    endforeach()
    set(_precision "")
    foreach(f ${${lib_name}_PRECISION_SOURCES})
        list(APPEND _precision "${_dir}/${f}")
    endforeach()

    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        fortran_module_layout(${_common_target})
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    fortran_module_layout(${_precision_target})
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()

    # Content-driven archive filename. The CMake target keeps LIB_PREFIX
    # (so exported symbols and inter-target references are unchanged), but
    # the emitted ``.a`` is renamed to reflect what it actually holds:
    # E/Y families → ey, Q/X → qx, M/W → mw (see ${lib_name}_ARCHIVE_PREFIX
    # in the staged manifest). fortran_install_library reads this property.
    set_target_properties(${_precision_target} PROPERTIES
        MIGRATED_OUTPUT_BASE "${${lib_name}_ARCHIVE_PREFIX}${lib_name}")

    # The migrator can lengthen tokens (e.g. MPI_DOUBLE_COMPLEX →
    # MPI_C_LONG_DOUBLE_COMPLEX for kind10), pushing fixed-form ``.F``
    # lines past column 72 and free-form ``.f90`` lines past column 132.
    # Disable both line-length limits so the post-migration source
    # compiles regardless of growth. (Free-form is latent today — the
    # widest migrated ``.f90`` sits at ~113 cols — but the same token
    # growth that overflows fixed form can overflow free form, so guard
    # both rather than wait for it to bite.)
    if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran"
       OR FORTRAN_COMPILER_FAMILY STREQUAL "flang")
        target_compile_options(${_precision_target} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>
            $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
        if(TARGET ${_common_target})
            target_compile_options(${_common_target} PRIVATE
                $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>
                $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
        endif()
    elseif(FORTRAN_COMPILER_FAMILY STREQUAL "intel")
        target_compile_options(${_precision_target} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-extend-source>)
        if(TARGET ${_common_target})
            target_compile_options(${_common_target} PRIVATE
                $<$<COMPILE_LANGUAGE:Fortran>:-extend-source>)
        endif()
    endif()

    # Link multifloats helpers for Fortran libraries. multifloatsf
    # comes from a fetched sub-project that owns its own install/export
    # set; wrap the link in ``$<BUILD_INTERFACE:>`` so the dependency
    # is visible during the build (tests, locally-built downstreams)
    # but isn't propagated into our install(EXPORT) — otherwise CMake
    # rejects the generate step because multifloatsf is "not in any
    # export set" known to us.
    if(NEEDS_MULTIFLOATS)
        target_link_libraries(${_precision_target} PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
        if(TARGET la_constants_mw)
            target_link_libraries(${_precision_target} PUBLIC la_constants_mw)
        endif()
        if(TARGET la_xisnan_mw)
            target_link_libraries(${_precision_target} PUBLIC la_xisnan_mw)
        endif()
        # The common archive can also hold family-independent sources that
        # the multifloats migrator rewrites to ``USE multifloats`` (e.g.
        # scalapack's PCHK1MAT/PCHK2MAT in pchkxmat.f rebind the generic
        # mod/min/operators onto TYPE(real64x2)). Those objects need the
        # multifloats module on their include path and a build-order
        # dependency on it. The precision target links multifloatsf above,
        # but that does not flow down to _common (precision links _common,
        # not the reverse), so link it here too. BUILD_INTERFACE keeps it
        # out of install(EXPORT), exactly as for the precision target.
        if(TARGET ${_common_target})
            target_link_libraries(${_common_target} PUBLIC
                $<BUILD_INTERFACE:multifloatsf>)
        endif()
    endif()
endfunction()


function(add_migrated_c_library lib_name)
    set(_dir ${CMAKE_CURRENT_SOURCE_DIR}/${lib_name})
    if(NOT EXISTS "${_dir}/manifest.cmake")
        return()
    endif()
    include(${_dir}/manifest.cmake)

    # Prefix source paths with library directory
    set(_common "")
    foreach(f ${${lib_name}_COMMON_SOURCES})
        list(APPEND _common "${_dir}/${f}")
    endforeach()
    set(_precision "")
    foreach(f ${${lib_name}_PRECISION_SOURCES})
        list(APPEND _precision "${_dir}/${f}")
    endforeach()
    set(_dual "")
    foreach(f ${${lib_name}_DUAL_INTERFACE_SOURCES})
        list(APPEND _dual "${_dir}/${f}")
    endforeach()

    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        target_include_directories(${_common_target} PUBLIC $<BUILD_INTERFACE:${_dir}/src>)
        if(MPI_C_FOUND)
            # Link MPI privately (needed for compiling the C sources) and
            # propagate it via INTERFACE link-only, so downstream Fortran
            # targets inherit the link dependency without MPI's C-only
            # compile flags (e.g. mpich's ``-ffat-lto-objects`` is rejected
            # by flang-new-20).
            target_link_libraries(${_common_target} PRIVATE MPI::MPI_C)
            target_link_libraries(${_common_target} INTERFACE $<LINK_ONLY:MPI::MPI_C>)
        endif()
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    target_include_directories(${_precision_target} PUBLIC $<BUILD_INTERFACE:${_dir}/src>)
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()
    if(MPI_C_FOUND)
        target_link_libraries(${_precision_target} PRIVATE MPI::MPI_C)
        target_link_libraries(${_precision_target} INTERFACE $<LINK_ONLY:MPI::MPI_C>)
    endif()

    # Content-driven archive filename (see add_migrated_fortran_library).
    set_target_properties(${_precision_target} PROPERTIES
        MIGRATED_OUTPUT_BASE "${${lib_name}_ARCHIVE_PREFIX}${lib_name}")

    # Dual-interface compile: files listed in
    # ``${lib_name}_DUAL_INTERFACE_SOURCES`` are already compiled once
    # above (they live in _precision / _common, since every such file
    # is also one of BLACS's regular sources) and expose the
    # Fortran-callable entry points by default. Compile them a second
    # time as an OBJECT library with ``-DCallFromC`` so the C-callable
    # entry points (``Cblacs_*``, ``Cdgesd2d`` etc.) are also present
    # in the final static archive.
    #
    # The C-callable variant must land in the SAME archive as its
    # Fortran-side counterpart, so a family-independent file (integer
    # ``igamx2d_`` or the agnostic driver) contributes its ``Cigamx2d`` /
    # ``Cblacs_*`` object to ``_common``, and only the e/y files
    # contribute ``Cegamx2d`` / ``Cygamx2d`` to the prefixed archive.
    # Otherwise the prefixed archive would re-acquire the integer /
    # agnostic C symbols and stop being strictly single-family. Partition
    # ``_dual`` by which regular-source list each file belongs to. Empty
    # partitions are a no-op.
    set(_dual_common "")
    set(_dual_precision "")
    foreach(f ${_dual})
        if("${f}" IN_LIST _common)
            list(APPEND _dual_common "${f}")
        else()
            list(APPEND _dual_precision "${f}")
        endif()
    endforeach()
    foreach(_slot common precision)
        if(_slot STREQUAL common)
            set(_dsrc "${_dual_common}")
            set(_dhost "${_common_target}")
        else()
            set(_dsrc "${_dual_precision}")
            set(_dhost "${_precision_target}")
        endif()
        if(NOT _dsrc OR NOT TARGET ${_dhost})
            continue()
        endif()
        set(_dual_target "${_dhost}_c_api")
        add_library(${_dual_target} OBJECT ${_dsrc})
        target_compile_definitions(${_dual_target} PRIVATE CallFromC)
        target_include_directories(${_dual_target} PRIVATE $<BUILD_INTERFACE:${_dir}/src>)
        if(MPI_C_FOUND)
            target_link_libraries(${_dual_target} PRIVATE MPI::MPI_C)
        endif()
        if(C_AS_CXX)
            set_source_files_properties(${_dsrc} PROPERTIES LANGUAGE CXX)
            target_compile_features(${_dual_target} PRIVATE cxx_std_17)
            target_compile_options(${_dual_target} PRIVATE -fpermissive -Wno-write-strings)
            # Skip the C++ MPI bindings (mpicxx.h). They define template
            # classes that conflict with the ``extern "C"`` wrap the
            # c_migrator post-pass injects around .c bodies. The C
            # bindings (defined in mpi.h itself) are all we need.
            target_compile_definitions(${_dual_target} PRIVATE
                MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
            # Multifloats: ``Bdef.h`` pulls in ``multifloats_bridge.h`` under
            # multifloats mode. That header is exposed through the
            # ``multifloats_mpi`` target (INTERFACE or STATIC depending on
            # MPI availability); propagate its include directories so the
            # dual-interface compile finds it just like the main target.
            if(TARGET multifloats_mpi)
                target_link_libraries(${_dual_target} PRIVATE multifloats_mpi)
            endif()
        endif()
        target_sources(${_dhost} PRIVATE $<TARGET_OBJECTS:${_dual_target}>)
    endforeach()

    # Multifloats: compile C as C++ for operator overloading support.
    # ``copy_files``-staged Fortran helpers (e.g. pblas/pilaenv.f) ride in
    # COMMON_SOURCES but must retain their native language; filter them
    # out of the language override.
    #
    # (ScaLAPACK's C-side REDIST sources are handled by the analogous
    # C_AS_CXX block in the scalapack_c fold below, not here — they are
    # never built through this function.)
    if(C_AS_CXX)
        set(_all_sources "")
        foreach(_s ${_common} ${_precision})
            get_filename_component(_ext "${_s}" EXT)
            string(TOLOWER "${_ext}" _ext)
            if(NOT _ext STREQUAL ".f" AND NOT _ext STREQUAL ".f90"
                AND NOT _ext STREQUAL ".f95" AND NOT _ext STREQUAL ".f03"
                AND NOT _ext STREQUAL ".for")
                list(APPEND _all_sources "${_s}")
            endif()
        endforeach()
        set_source_files_properties(${_all_sources} PROPERTIES LANGUAGE CXX)
        if(TARGET ${_common_target})
            target_compile_features(${_common_target} PRIVATE cxx_std_17)
            target_compile_options(${_common_target} PRIVATE -fpermissive -Wno-write-strings)
            target_compile_definitions(${_common_target} PRIVATE
                MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
        endif()
        target_compile_features(${_precision_target} PRIVATE cxx_std_17)
        target_compile_options(${_precision_target} PRIVATE -fpermissive -Wno-write-strings)
        target_compile_definitions(${_precision_target} PRIVATE
            MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
        # Link multifloats + bridge for C libraries in multifloats mode
        if(TARGET multifloats)
            target_link_libraries(${_precision_target} PUBLIC multifloats)
        endif()
        if(TARGET multifloats_mpi)
            target_link_libraries(${_precision_target} PUBLIC multifloats_mpi)
            if(TARGET ${_common_target})
                target_link_libraries(${_common_target} PUBLIC multifloats_mpi)
            endif()
        endif()
    endif()
endfunction()

# ── Libraries in dependency order ────────────────────────────────────

# 1. BLAS (no deps).
#
# Standard-precision archive ``blas`` is built from the vendored Netlib
# BLAS sources (staged into ``_refblas_src/`` by ``migrator stage``).
# The migrated archive ``${LIB_PREFIX}blas`` carries only Q/X/E/Y/M/W
# clones and PUBLIC-links the standard archive so downstreams resolve
# the S/D/C/Z entry points (DGEMM_, ZGEMM_, ...) through one link line.
# Cross-precision routines (dsdot, sdsdot) are excluded — the migrator
# also skips them (see recipes/blas.yaml).
add_standard_fortran_library(blas
    ${CMAKE_CURRENT_SOURCE_DIR}/_refblas_src
    EXCLUDE_REGEX "/(dsdot|sdsdot)\\.f9?0?$")
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(blas)
    if(TARGET ${LIB_PREFIX}blas AND TARGET blas)
        target_link_libraries(${LIB_PREFIX}blas PUBLIC blas)
    endif()
endif()

# 1.5. XBLAS (depends on BLAS).
#
# Netlib's Extra-Precise BLAS (1.0.248) is a C library that exposes
# Level 1/2/3 BLAS-like entry points with an additional ``prec``
# parameter selecting the internal accumulation precision
# (single / double / extra). The migrated archive ``${LIB_PREFIX}xblas``
# carries Q/X/E/Y clones (BLAS_qgemv_x, BLAS_xgemv_x, ...) of every
# uniform-precision _x routine; mixed-input variants
# (BLAS_dgemv_d_s_x, ...) are skipped by the recipe — they're not
# expressible in a single extended-precision target and unreachable
# from LAPACK.
#
# There is no standard-precision archive for XBLAS in this build:
# nothing in the project links the regular C XBLAS today (LAPACK's
# *_extended.f routines that *would* call into it are deferred to a
# follow-up that adds a -DBUILD_XBLAS gate). We can revive a
# standard archive when that need lands.
if(NOT BASELINE_BUILD)
    add_migrated_c_library(xblas)
    if(TARGET ${LIB_PREFIX}xblas)
        # f2c-bridge.h defaults to no Fortran-symbol mangling. gfortran
        # emits lowercase + trailing underscore (the "Add_" convention),
        # so define CONFIG_FC_UNDERSCORE to make the bridge produce
        # ``blas_qgemv_x_`` etc. — matching what Fortran callers expect.
        target_compile_definitions(${LIB_PREFIX}xblas PRIVATE
            CONFIG_FC_UNDERSCORE)
        if(TARGET xblas_common)
            target_compile_definitions(xblas_common PRIVATE
                CONFIG_FC_UNDERSCORE)
        endif()
        if(TARGET ${LIB_PREFIX}blas)
            target_link_libraries(${LIB_PREFIX}xblas PUBLIC ${LIB_PREFIX}blas)
        endif()
    endif()
endif()

# 2. BLACS (no library deps, needs MPI).
#
# Standard archive ``blacs`` is built from the vendored upstream BLACS
# C sources (staged into ``_blacs_src/``). ``-DAdd_`` selects the
# Fortran-symbol-mangling convention gfortran emits (lowercase + trailing
# underscore), matching what the migrated archive expects. The build
# only stages the Fortran-callable entry points; the dual C-callable
# entry points (Cblacs_*, Cdgesd2d, ...) needed by the migrated
# scalapack C clones (folded into the migrated scalapack archive) are
# produced by the migrated archive's DUAL_INTERFACE_SOURCES path — the
# std archive doesn't replicate them.
add_standard_c_library(blacs
    ${CMAKE_CURRENT_SOURCE_DIR}/_blacs_src
    COMPILE_DEFINITIONS Add_)
# BLACS sources gate their entry-point signature on
# ``INTFACE == C_CALL`` / ``CallFromC`` so each .c file produces two
# distinct ABI surfaces: the Fortran-callable ``blacs_*_`` set
# (default compile) and the C-callable ``Cblacs_*`` set
# (-DCallFromC). The migrated archive's add_migrated_c_library
# handles this via DUAL_INTERFACE_SOURCES; mirror it for the std
# archive by compiling the same sources a second time with
# CallFromC defined and folding the resulting objects into ``blacs``.
if(TARGET blacs)
    file(GLOB _blacs_dual_srcs CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/_blacs_src/*.c)
    add_library(blacs_c_api OBJECT ${_blacs_dual_srcs})
    target_compile_definitions(blacs_c_api PRIVATE Add_ CallFromC)
    target_include_directories(blacs_c_api PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_blacs_src>)
    if(MPI_C_FOUND)
        target_link_libraries(blacs_c_api PRIVATE MPI::MPI_C)
    endif()
    target_sources(blacs PRIVATE $<TARGET_OBJECTS:blacs_c_api>)
endif()
if(NOT BASELINE_BUILD)
    add_migrated_c_library(blacs)
    if(TARGET ${LIB_PREFIX}blacs AND NEEDS_MULTIFLOATS AND TARGET multifloats_mpi)
        target_link_libraries(${LIB_PREFIX}blacs PUBLIC multifloats_mpi)
    endif()
    if(TARGET ${LIB_PREFIX}blacs AND NEEDS_QUAD_MPI AND TARGET quad_mpi)
        target_link_libraries(${LIB_PREFIX}blacs PUBLIC quad_mpi)
    endif()
    if(TARGET ${LIB_PREFIX}blacs AND TARGET blacs)
        target_link_libraries(${LIB_PREFIX}blacs PUBLIC blacs)
    endif()
endif()

# 3. LAPACK (depends on BLAS).
#
# Standard archive ``lapack`` is built from the vendored Netlib LAPACK
# sources (staged into ``_reflapack_src/``). Mixed-precision bridge
# routines (dsgesv, zcgesv, dsposv, zcposv, dsgels, zcgels, dlag2s,
# slag2d, zlag2c, clag2z, dlat2s, zlat2c) are excluded — they don't
# promote under -freal-8-real-16, but the whole point of the standard
# archive is to NOT promote, so these are simply left out for parity
# with the migrated archive's skip list. Extended-precision-only
# helper modules (la_constants_{ey,qx,mw} / la_xisnan_{ey,qx,mw}) are
# excluded — they only make sense in the migrated archive since they
# reference extended / quad / multifloats kinds.
add_standard_fortran_library(lapack
    ${CMAKE_CURRENT_SOURCE_DIR}/_reflapack_src
    EXCLUDE_REGEX "/(dsgesv|zcgesv|dsposv|zcposv|dsgels|zcgels|dlag2s|slag2d|zlag2c|clag2z|dlat2s|zlat2c|la_constants_ey|la_constants_qx|la_constants_mw|la_xisnan_ey|la_xisnan_qx|la_xisnan_mw)\\.([fF]|[fF]90)$")
if(TARGET lapack AND TARGET blas)
    target_link_libraries(lapack PUBLIC blas)
endif()
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(lapack)
    if(TARGET ${LIB_PREFIX}lapack AND TARGET ${LIB_PREFIX}blas)
        target_link_libraries(${LIB_PREFIX}lapack PUBLIC ${LIB_PREFIX}blas)
    endif()
    if(TARGET ${LIB_PREFIX}lapack AND TARGET lapack)
        target_link_libraries(${LIB_PREFIX}lapack PUBLIC lapack)
    endif()
endif()
# The xx-family iterative-refinement helpers (qla_gerfsx_extended,
# xla_gerfsx_extended, ...) live inside ${LIB_PREFIX}lapack and call
# XBLAS extra-precise Level-2 entry points (BLAS_QGEMV_X_, ...). With
# static archives, linking ${LIB_PREFIX}lapack into a test executable
# pulls those *_extended objects via undefined-symbol resolution as
# soon as any xx driver (qgesvxx / xgesvxx / qgbsvxx / ...) is reached
# transitively, so the XBLAS symbols must be resolvable on the same
# link line. Make the xblas dependency public so every consumer of
# the migrated lapack archive picks it up automatically.
if(NOT BASELINE_BUILD)
    if(TARGET ${LIB_PREFIX}lapack AND TARGET ${LIB_PREFIX}xblas)
        target_link_libraries(${LIB_PREFIX}lapack PUBLIC ${LIB_PREFIX}xblas)
    endif()
    # The per-target precision helpers (la_constants_{ey,qx},
    # la_xisnan_{ey,qx}) now ride in the PREFIXED ${LIB_PREFIX}lapack
    # archive (routed there by staging.py), not the shared lapack_common.
    # They ``USE`` the standard la_constants / la_xisnan modules supplied
    # by the std ``lapack`` archive; the ${LIB_PREFIX}lapack -> lapack
    # PUBLIC link above already resolves ``la_xisnan.mod`` at compile
    # time. The lapack_common -> lapack link below is now defensive — it
    # only fires if some other precision-independent member keeps a
    # lapack_common archive alive.
    if(TARGET lapack_common AND TARGET lapack)
        target_link_libraries(lapack_common PUBLIC lapack)
    endif()
endif()

# 4. PTZBLAS (depends on BLAS).
#
# Standard archive ``ptzblas`` is built from PBLAS/PTZBLAS/. Fortran
# kernels reference ZZDOTC/ZZDOTU (defined in scalapack/TOOLS/) — those
# resolve through std ``scalapack`` at final link time, mirroring the
# upstream layering.
add_standard_fortran_library(ptzblas
    ${CMAKE_CURRENT_SOURCE_DIR}/_ptzblas_src)
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(ptzblas)
    if(TARGET ${LIB_PREFIX}ptzblas AND TARGET ${LIB_PREFIX}blas)
        target_link_libraries(${LIB_PREFIX}ptzblas PUBLIC ${LIB_PREFIX}blas)
    endif()
    if(TARGET ${LIB_PREFIX}ptzblas AND TARGET ptzblas)
        target_link_libraries(${LIB_PREFIX}ptzblas PUBLIC ptzblas)
    endif()
endif()

# 4b. ScaLAPACK tools (NUMROC / ICEIL / ILCM) — FOLDED INTO scalapack_common.
#
# These three precision-independent integer helpers used to form a
# standalone ``${LIB_PREFIX}scalapack_tools`` archive, built once per
# family with a prefix. That duplicated three symbols MKL also provides
# and made them a mandatory dep of ``${LIB_PREFIX}pbblas``, so under
# extended-MKL they multiply-defined against MKL's numroc_/iceil_/ilcm_
# — the exact overlap Task 29 removes. They are now migrated by
# scalapack.yaml into the family-independent ``scalapack_common`` (which
# the consumer drops when MKL supplies that layer, same profile as
# ``blacs_common``). PBBLAS reaches them by PUBLIC-linking
# ``scalapack_common`` once it exists (wired after section 8, below).
# The standard ``scalapack`` archive still absorbs TOOLS/*.f directly
# via _scalapack_tools_src (see section 7) — unchanged.

# 5. PBBLAS (depends on BLAS, BLACS; scalapack_common wired post-section-8).
#
# Standard archive ``pbblas`` is built from PBLAS/PBBLAS/. Pure Fortran;
# no extra dependencies beyond BLAS+BLACS at link time.
add_standard_fortran_library(pbblas
    ${CMAKE_CURRENT_SOURCE_DIR}/_pbblas_src)
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(pbblas)
    if(TARGET ${LIB_PREFIX}pbblas)
        foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}pbblas PUBLIC ${_dep})
            endif()
        endforeach()
        # NUMROC/ICEIL/ILCM now live in scalapack_common (folded in from
        # the old scalapack_tools archive). scalapack_common is created
        # later, in section 7's add_migrated_fortran_library(scalapack),
        # so the PUBLIC link is wired after section 8 where the target
        # exists — see the "PBBLAS → scalapack_common" block below.
    endif()
    if(TARGET ${LIB_PREFIX}pbblas AND TARGET pbblas)
        target_link_libraries(${LIB_PREFIX}pbblas PUBLIC pbblas)
    endif()
endif()

# 6. PBLAS (depends on BLAS, BLACS, PTZBLAS, PBBLAS; C library).
#
# Standard archive ``pblas`` is built from PBLAS/SRC/*.c plus
# PBLAS/SRC/PTOOLS/*.c (the helper directory the migrated archive
# pulls in via ``extra_c_dirs``). pilaenv.f (a precision-independent
# Fortran block-size query) is added too. ``-DAdd_`` selects the
# Fortran symbol-mangling convention.
add_standard_c_library(pblas
    ${CMAKE_CURRENT_SOURCE_DIR}/_pblas_src
    COMPILE_DEFINITIONS Add_
    INCLUDE_DIRS
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_pblas_src>)
if(TARGET pblas)
    # PTOOLS sources live under _pblas_src/PTOOLS/ (matching upstream's
    # PBLAS/SRC/PTOOLS layout) so ``#include "../pblas.h"`` resolves
    # relative to the file's own directory back into _pblas_src/.
    file(GLOB _pblas_ptools_srcs CONFIGURE_DEPENDS
        ${CMAKE_CURRENT_SOURCE_DIR}/_pblas_src/PTOOLS/*.c)
    if(_pblas_ptools_srcs)
        target_sources(pblas PRIVATE ${_pblas_ptools_srcs})
    endif()
    set(_pblas_pilaenv ${CMAKE_CURRENT_SOURCE_DIR}/_pblas_src/pilaenv.f)
    if(EXISTS ${_pblas_pilaenv})
        target_sources(pblas PRIVATE ${_pblas_pilaenv})
        # pilaenv.f is fixed-form Fortran; CMake auto-detects from .f.
    endif()
endif()
if(NOT BASELINE_BUILD)
    add_migrated_c_library(pblas)
    if(TARGET ${LIB_PREFIX}pblas)
        foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs
                     ${LIB_PREFIX}ptzblas ${LIB_PREFIX}pbblas)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}pblas PUBLIC ${_dep})
            endif()
        endforeach()
    endif()
    if(TARGET ${LIB_PREFIX}pblas AND TARGET pblas)
        target_link_libraries(${LIB_PREFIX}pblas PUBLIC pblas)
    endif()
endif()

# 7. ScaLAPACK (depends on LAPACK, BLACS, PBLAS).
#
# Standard archive ``scalapack`` is built from SRC/*.f (Fortran entry
# points) plus TOOLS/*.f / *.f90 (descriptor/grid/communication
# helpers — descinit_, numroc_, iceil_, the pd*/pz* element accessors,
# zzdotc_/zzdotu_, ...). Mixed-precision iterative drivers are
# excluded for parity with the migrated recipe — they're not
# expressible in extended precision and the std archive has no need
# to ship them either.
add_standard_fortran_library(scalapack
    ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src
    EXCLUDE_REGEX "/(zcgesv|dsgesv|zcposv|dsposv|zcgels|dsgels|dlag2s|slag2d|clag2z|zlag2c)\\.[fF]9?0?$")
if(TARGET scalapack)
    file(GLOB _scalapack_tools_f
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src/*.f)
    file(GLOB _scalapack_tools_f90
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src/*.f90)
    set(_scalapack_tools ${_scalapack_tools_f} ${_scalapack_tools_f90})
    if(_scalapack_tools)
        target_sources(scalapack PRIVATE ${_scalapack_tools})
    endif()
endif()
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(scalapack)
    if(TARGET ${LIB_PREFIX}scalapack)
        # scalapack_common (which now carries NUMROC/ICEIL/ILCM) is linked
        # automatically by add_migrated_fortran_library(scalapack) above.
        foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}scalapack PUBLIC ${_dep})
            endif()
        endforeach()
    endif()
    if(TARGET ${LIB_PREFIX}scalapack AND TARGET scalapack)
        target_link_libraries(${LIB_PREFIX}scalapack PUBLIC scalapack)
    endif()
endif()

# 8. ScaLAPACK C entry points — FOLDED INTO the scalapack archive.
#
# The C entry points under SRC/*.c (clamov, dlamov, pdlaiect, pcrot,
# ...) plus REDIST/SRC/*.c (pdgemr2d_, pztrmr2d_, ...) and TOOLS/*.c
# (sl_gridreshape, reshape) used to form a standalone ``scalapack_c``
# archive that PUBLIC-linked ``scalapack``, while ``scalapack``'s own
# Fortran routines call back into these C helpers (qlamov_, pqgemr2d_).
# That mutual dependency was a deliberately-broken cycle: wiring the
# reverse link put qscalapack_c into qscalapack's
# INTERFACE_LINK_LIBRARIES, and the resulting cycle broke
# ``find_package`` (qscalapack's Targets file failed the
# ``_cmake_expected_target`` check before qscalapack_c loaded), so
# consumers had to link BOTH archives explicitly.
#
# Folding the C objects directly into the scalapack archive makes the
# reference INTERNAL to one archive: ``ld`` resolves intra-archive
# back-references via the archive symbol index regardless of member
# order (no ``--start-group`` needed), so the cycle disappears and a
# single ``eplinalg::scalapack`` carries both halves. Built as OBJECT
# libraries to keep the C (and, under multifloats, C-as-C++) compile
# settings off the Fortran sources. ``-DAdd_`` selects the lowercase +
# trailing-underscore Fortran symbol mangling that gfortran emits
# (upstream SLmake.inc.example's convention).
if(TARGET scalapack)
    file(GLOB _scalapack_c_std
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src/*.c)
    file(GLOB _scalapack_c_redist
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_redist_src/*.c)
    file(GLOB _scalapack_c_tools
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src/*.c)
    set(_scalapack_c_all
        ${_scalapack_c_std} ${_scalapack_c_redist} ${_scalapack_c_tools})
    if(_scalapack_c_all)
        add_library(scalapack_c_std_obj OBJECT ${_scalapack_c_all})
        target_compile_definitions(scalapack_c_std_obj PRIVATE Add_)
        target_include_directories(scalapack_c_std_obj PRIVATE
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_redist_src>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src>)
        if(MPI_C_FOUND)
            target_link_libraries(scalapack_c_std_obj PRIVATE MPI::MPI_C)
        endif()
        target_sources(scalapack PRIVATE $<TARGET_OBJECTS:scalapack_c_std_obj>)
    endif()
    # scalapack now carries C objects that call MPI (the standalone
    # scalapack_c archive used to pull it); propagate the link privately
    # for its own objects and interface-only for downstream Fortran.
    if(MPI_C_FOUND)
        target_link_libraries(scalapack PRIVATE MPI::MPI_C)
        target_link_libraries(scalapack INTERFACE $<LINK_ONLY:MPI::MPI_C>)
    endif()
endif()
if(NOT BASELINE_BUILD AND TARGET ${LIB_PREFIX}scalapack)
    # Migrated C clones (qlamov_, pqgemr2d_, ...) fold into the migrated
    # scalapack archive pair instead of a separate qscalapack_c. Read the
    # scalapack_c manifest for the precision/common C split and build each
    # as an OBJECT lib; the precision objects fold into
    # ${LIB_PREFIX}scalapack and the arithmetic-independent common object
    # (pbchkvect.c) into scalapack_common.
    set(_scalapack_c_dir ${CMAKE_CURRENT_SOURCE_DIR}/scalapack_c)
    if(EXISTS "${_scalapack_c_dir}/manifest.cmake")
        include(${_scalapack_c_dir}/manifest.cmake)
        set(_scalapack_c_prec_src "")
        foreach(f ${scalapack_c_PRECISION_SOURCES})
            list(APPEND _scalapack_c_prec_src "${_scalapack_c_dir}/${f}")
        endforeach()
        set(_scalapack_c_common_src "")
        foreach(f ${scalapack_c_COMMON_SOURCES})
            list(APPEND _scalapack_c_common_src "${_scalapack_c_dir}/${f}")
        endforeach()
        if(_scalapack_c_prec_src)
            add_library(scalapack_c_prec_obj OBJECT ${_scalapack_c_prec_src})
        endif()
        if(_scalapack_c_common_src)
            add_library(scalapack_c_common_obj OBJECT ${_scalapack_c_common_src})
        endif()
        foreach(_obj scalapack_c_prec_obj scalapack_c_common_obj)
            if(TARGET ${_obj})
                target_compile_definitions(${_obj} PRIVATE Add_)
                target_include_directories(${_obj} PRIVATE
                    $<BUILD_INTERFACE:${_scalapack_c_dir}/src>)
                if(MPI_C_FOUND)
                    target_link_libraries(${_obj} PRIVATE MPI::MPI_C)
                endif()
            endif()
        endforeach()
        if(C_AS_CXX)
            # Multifloats: compile as C++ (same treatment as
            # add_migrated_c_library). The REDIST K&R-style definitions the
            # migrator leaves unconverted are tolerated by ``-fpermissive``;
            # these routines do only block-cyclic copies and sizeof pointer
            # arithmetic — never per-element multifloats math — so they need
            # no operator overloading, but still ride the CXX override for a
            # uniform multifloats C ABI.
            set_source_files_properties(
                ${_scalapack_c_prec_src} ${_scalapack_c_common_src}
                PROPERTIES LANGUAGE CXX)
            foreach(_obj scalapack_c_prec_obj scalapack_c_common_obj)
                if(TARGET ${_obj})
                    target_compile_features(${_obj} PRIVATE cxx_std_17)
                    target_compile_options(${_obj} PRIVATE
                        -fpermissive -Wno-write-strings)
                    target_compile_definitions(${_obj} PRIVATE
                        MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
                    if(TARGET multifloats_mpi)
                        target_link_libraries(${_obj} PRIVATE multifloats_mpi)
                    endif()
                endif()
            endforeach()
        endif()
        if(TARGET scalapack_c_prec_obj)
            target_sources(${LIB_PREFIX}scalapack
                PRIVATE $<TARGET_OBJECTS:scalapack_c_prec_obj>)
        endif()
        if(TARGET scalapack_c_common_obj)
            if(TARGET scalapack_common)
                target_sources(scalapack_common
                    PRIVATE $<TARGET_OBJECTS:scalapack_c_common_obj>)
            else()
                target_sources(${LIB_PREFIX}scalapack
                    PRIVATE $<TARGET_OBJECTS:scalapack_c_common_obj>)
            endif()
        endif()
        # The migrated scalapack archive now carries the C objects: expose
        # the C source dir (internal prototype headers) and propagate MPI
        # + the multifloats bridge the standalone archive used to.
        target_include_directories(${LIB_PREFIX}scalapack
            PUBLIC $<BUILD_INTERFACE:${_scalapack_c_dir}/src>)
        if(MPI_C_FOUND)
            target_link_libraries(${LIB_PREFIX}scalapack PRIVATE MPI::MPI_C)
            target_link_libraries(${LIB_PREFIX}scalapack
                INTERFACE $<LINK_ONLY:MPI::MPI_C>)
        endif()
        if(C_AS_CXX)
            if(TARGET multifloats)
                target_link_libraries(${LIB_PREFIX}scalapack PUBLIC multifloats)
            endif()
            if(TARGET multifloats_mpi)
                target_link_libraries(${LIB_PREFIX}scalapack PUBLIC multifloats_mpi)
            endif()
        endif()
    endif()
endif()

# PBBLAS → scalapack_common (NUMROC / ICEIL / ILCM).
#
# pb*tran / pb*trnv / pb*trsrt / pb*trget / pb*trst1 in PBBLAS call
# NUMROC / ICEIL / ILCM. Those three precision-independent helpers were
# folded from the retired ${LIB_PREFIX}scalapack_tools archive into
# scalapack_common (section 4b). scalapack_common is created in
# section 7 (after PBBLAS, section 5), so the PUBLIC link is wired here
# where the target is guaranteed to exist. tests/pbblas links only
# ${LIB_PREFIX}pbblas — it now transitively pulls scalapack_common for
# these symbols instead of the old tools archive, without dragging in
# the full ${LIB_PREFIX}scalapack.
if(NOT BASELINE_BUILD AND TARGET ${LIB_PREFIX}pbblas AND TARGET scalapack_common)
    target_link_libraries(${LIB_PREFIX}pbblas PUBLIC scalapack_common)
endif()

# 9. MUMPS (depends on LAPACK, BLAS, ScaLAPACK; needs MPI Fortran).
#
# MUMPS is a sparse direct solver (Multifrontal Massively Parallel). The
# migrated archive ``${LIB_PREFIX}mumps`` carries the Q/X clones of
# DMUMPS / ZMUMPS plus their auxiliary subroutines.
#
# No standard-precision archive is built today: nothing in this project
# links into the upstream-precision DMUMPS_/ZMUMPS_ entry points, and
# building one would require wiring the C side (mumps_c.c, mumps_common.c,
# mumps_pord.c, ...), which the recipe currently skips. Add later if a
# real consumer surfaces.
#
# MPI is linked PUBLIC because the migrated routines invoke MPI
# primitives (MPI_BCAST / MPI_REDUCE / ...). Downstream callers pass
# ``id%COMM = MPI_COMM_WORLD`` and let MUMPS handle communication
# internally; the link dependency must propagate so executables resolve
# the MPI runtime.
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(mumps)

    # ── PORD ordering library (ships in-tree with MUMPS) ────────────
    # PORD is a self-contained standard-C nested-dissection ordering
    # (no MPI, no external dependency); its 14 algorithm sources live
    # under MUMPS_5.8.2/PORD/lib and are staged verbatim into
    # _mumps_pord_src. Building ``libpord`` here and defining ``-Dpord``
    # on the MUMPS C runtime + Fortran analysis routines activates the
    # ICNTL(7)=4 ordering; without the define, mumps_pord.c compiles as
    # an inert stub (its body is guarded by ``#if defined(pord)``) and
    # the Fortran never references PORD, so the archive is unused dead
    # code. PORD permutes the INTEGER adjacency graph — it never touches
    # FP values — so this single build links into every migrated
    # arithmetic (and the genuine double MUMPS) alike.
    set(_mumps_pord_src ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_pord_src)
    set(_mumps_pord_inc ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_pord_include)
    set(MUMPS_HAVE_PORD FALSE)
    if(IS_DIRECTORY ${_mumps_pord_src} AND IS_DIRECTORY ${_mumps_pord_inc})
        file(GLOB _pord_c CONFIGURE_DEPENDS ${_mumps_pord_src}/*.c)
        if(_pord_c)
            add_library(pord STATIC ${_pord_c})
            # PORD symbols (SPACE_ordering, firstPostorder, …) are called
            # by plain C name from mumps_pord.c — no Fortran mangling, so
            # no ``Add_`` here. Its sources #include <space.h> etc. from
            # PORD/include (angle-bracket), hence the PUBLIC include dir,
            # which also reaches mumps_pord.c once the runtime links pord.
            target_include_directories(pord PUBLIC
                $<BUILD_INTERFACE:${_mumps_pord_inc}>)
            set(MUMPS_HAVE_PORD TRUE)
            message(STATUS "MUMPS: PORD ordering enabled (ICNTL(7)=4)")
        endif()
    endif()
    if(NOT MUMPS_HAVE_PORD)
        message(STATUS "MUMPS: PORD ordering unavailable "
                       "(no _mumps_pord_src) — ICNTL(7)=4 will be inactive")
    endif()

    # ── METIS ordering library (privately namespaced) ───────────────
    # METIS 5.1.0, vendored under external/metis-5.1.0 and staged into
    # _mumps_metis_{gklib,lib,include}. Every public API symbol was
    # renamed METIS_<X> → METIS_MUMPS_<X> (and internal libmetis__ →
    # libmetis_MUMPS_) so this private copy can never clash at link time
    # with a system METIS the final application might also pull in; the
    # matching MUMPS caller sites in _mumps_upstream_src were renamed too.
    # Building ``libmetis`` here and defining ``-Dmetis`` on the MUMPS C
    # runtime + Fortran analysis routines activates the ICNTL(7)=5
    # ordering; without the define mumps_metis*.c compile as inert stubs
    # (bodies guarded by ``#if defined(metis)``) and the Fortran never
    # references METIS. METIS permutes the INTEGER adjacency graph — it
    # never touches FP values — so a single 32-bit-idx build
    # (IDXTYPEWIDTH=32) links into every migrated arithmetic and the
    # genuine double MUMPS alike.
    set(_mumps_metis_gklib ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_metis_gklib)
    set(_mumps_metis_lib   ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_metis_lib)
    set(_mumps_metis_inc   ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_metis_include)
    set(MUMPS_HAVE_METIS FALSE)
    if(IS_DIRECTORY ${_mumps_metis_gklib} AND IS_DIRECTORY ${_mumps_metis_lib}
            AND IS_DIRECTORY ${_mumps_metis_inc})
        file(GLOB _metis_c CONFIGURE_DEPENDS
            ${_mumps_metis_gklib}/*.c ${_mumps_metis_lib}/*.c)
        if(_metis_c)
            add_library(metis STATIC ${_metis_c})
            # GKlib.h / metis.h / rename.h / defs.h … are angle-bracket
            # includes resolved from all three staged dirs; PUBLIC so the
            # public metis.h also reaches mumps_metis*.c once the runtime
            # links metis.
            target_include_directories(metis PUBLIC
                $<BUILD_INTERFACE:${_mumps_metis_inc}>
                $<BUILD_INTERFACE:${_mumps_metis_gklib}>
                $<BUILD_INTERFACE:${_mumps_metis_lib}>)
            # Width macros: 32-bit idx / real — the header default and the
            # width mumps_metis.c's ``IDXTYPEWIDTH == 32`` path expects.
            # NDEBUG* match the upstream (non-ASSERT) GKlib build; the
            # platform/feature macros mirror GKlib/GKlibSystem.cmake.
            target_compile_definitions(metis PRIVATE
                IDXTYPEWIDTH=32 REALTYPEWIDTH=32
                NDEBUG NDEBUG2 _FILE_OFFSET_BITS=64)
            if(UNIX AND NOT APPLE)
                target_compile_definitions(metis PRIVATE LINUX)
            endif()
            include(CheckIncludeFile)
            include(CheckFunctionExists)
            check_include_file(execinfo.h _METIS_HAVE_EXECINFO_H)
            if(_METIS_HAVE_EXECINFO_H)
                target_compile_definitions(metis PRIVATE HAVE_EXECINFO_H)
            endif()
            check_function_exists(getline _METIS_HAVE_GETLINE)
            if(_METIS_HAVE_GETLINE)
                target_compile_definitions(metis PRIVATE HAVE_GETLINE)
            endif()
            set_target_properties(metis PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON)
            # GKlib is warning-noisy by design (upstream builds it with a
            # long -Wno-* list); quiet it rather than drown the log.
            if(NOT MSVC)
                target_compile_options(metis PRIVATE -w -fno-strict-aliasing)
            endif()
            set(MUMPS_HAVE_METIS TRUE)
            message(STATUS "MUMPS: METIS ordering enabled (ICNTL(7)=5)")
        endif()
    endif()
    if(NOT MUMPS_HAVE_METIS)
        message(STATUS "MUMPS: METIS ordering unavailable "
                       "(no _mumps_metis_lib) — ICNTL(7)=5 will be inactive")
    endif()

    # ── Scotch ordering library (privately namespaced) ──────────────
    # Scotch 7.0.4, vendored under external/scotch-7.0.4 and staged into
    # _mumps_scotch_{libsrc,esmumps,include}. The archives are compiled
    # with -DSCOTCH_NAME_SUFFIX=_mumps so every public SCOTCH_* and
    # internal _SCOTCH* symbol carries a _mumps suffix — this private copy
    # can never clash at link time with a system Scotch the final
    # application might pull in. The bison/flex parser (parser_yy.c,
    # parser_ll.c, parser_ly.h) and the public headers (scotch.h,
    # scotchf.h) are pre-generated and vendored, so the eplinalg build
    # gains no bison/flex dependency. Building ``scotch`` (+ ``scotcherr``,
    # ``esmumps`` the MUMPS adapter) here and defining ``-Dscotch`` on the
    # MUMPS C runtime turns on the ICNTL(7)=3 sequential ordering; without
    # the define mumps_scotch*.c compile as inert stubs (bodies guarded by
    # ``#if defined(scotch)``) and the Fortran never references Scotch.
    # Scotch permutes the INTEGER adjacency graph — precision-agnostic, so
    # a single build serves every migrated arithmetic. PT-Scotch stays off
    # (sequential ICNTL(7)=3 only). The MUMPS caller sites reference bare
    # SCOTCH_* names; scotch_rename_mumps.h (force-included on the C
    # runtime and the test bridge) maps them onto the suffixed symbols.
    set(_mumps_scotch_libsrc ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_scotch_libsrc)
    set(_mumps_scotch_esmsrc ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_scotch_esmumps)
    set(_mumps_scotch_inc    ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_scotch_include)
    set(MUMPS_HAVE_SCOTCH FALSE)
    if(IS_DIRECTORY ${_mumps_scotch_libsrc} AND IS_DIRECTORY ${_mumps_scotch_esmsrc}
            AND IS_DIRECTORY ${_mumps_scotch_inc}
            AND EXISTS ${_mumps_scotch_inc}/scotch_sources.cmake)
        # Authoritative compile lists (basenames): PT-Scotch/distributed
        # .c live in the tree for textual #include but are not compiled.
        include(${_mumps_scotch_inc}/scotch_sources.cmake)
        set(_scotch_lib_c "")
        foreach(_f ${SCOTCH_LIBSCOTCH_SOURCES})
            list(APPEND _scotch_lib_c ${_mumps_scotch_libsrc}/${_f})
        endforeach()
        set(_scotch_err_c "")
        foreach(_f ${SCOTCH_SCOTCHERR_SOURCES})
            list(APPEND _scotch_err_c ${_mumps_scotch_libsrc}/${_f})
        endforeach()
        set(_scotch_esm_c "")
        foreach(_f ${SCOTCH_ESMUMPS_SOURCES})
            list(APPEND _scotch_esm_c ${_mumps_scotch_esmsrc}/${_f})
        endforeach()
        if(_scotch_lib_c AND _scotch_esm_c)
            # Exact define set the upstream Scotch CMake applies, plus the
            # namespacing suffix. -Drestrict=__restrict quiets the C99
            # keyword; COMMON_RANDOM_FIXED_SEED makes orderings reproducible.
            set(_scotch_defs
                COMMON_RANDOM_FIXED_SEED
                SCOTCH_VERSION_NUM=7 SCOTCH_RELEASE_NUM=0 SCOTCH_PATCHLEVEL_NUM=4
                SCOTCH_RENAME SCOTCH_NAME_SUFFIX=_mumps
                restrict=__restrict NDEBUG)

            add_library(scotch STATIC ${_scotch_lib_c})
            target_include_directories(scotch PUBLIC
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>
                $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
            target_compile_definitions(scotch PRIVATE ${_scotch_defs})
            set_target_properties(scotch PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON)

            add_library(scotcherr STATIC ${_scotch_err_c})
            target_include_directories(scotcherr PUBLIC
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>
                $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
            target_compile_definitions(scotcherr PRIVATE ${_scotch_defs})
            set_target_properties(scotcherr PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON)

            add_library(esmumps STATIC ${_scotch_esm_c})
            target_include_directories(esmumps PUBLIC
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
            target_include_directories(esmumps PRIVATE
                $<BUILD_INTERFACE:${_mumps_scotch_esmsrc}>
                $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
            target_compile_definitions(esmumps PRIVATE ${_scotch_defs})
            set_target_properties(esmumps PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON)
            # esmumps.c / order_scotch_graph.c / … include the *public*
            # scotch.h (bare prototypes) and call bare SCOTCH_* — its own
            # esmumps/module.h does NOT carry the public rename list (only
            # libscotch/module.h does). Force-include scotch_rename_mumps.h
            # so those bare calls bind to the _mumps-suffixed archive
            # symbols, exactly as the MUMPS C callers do.
            target_compile_options(esmumps PRIVATE
                "SHELL:-include ${_mumps_scotch_inc}/scotch_rename_mumps.h")
            target_link_libraries(esmumps PUBLIC scotch)

            # Scotch is warning-noisy; quiet it rather than drown the log.
            # It also leans on a "based array" idiom (`tab = ptr - baseval`,
            # a pointer legally *before* the allocation, re-based on every
            # access). Combined with the distro-default _FORTIFY_SOURCE and
            # __builtin_dynamic_object_size at -O2/-O3, that negative-offset
            # base defeats the compiler's object-size analysis and yields a
            # FALSE "buffer overflow detected" abort in symbolFax's in-bounds
            # memset (esmumps/symbol_fax.c:281). Upstream Scotch never sees
            # this because its own Makefiles do not inject fortify; match that
            # by turning it off for the vendored sources only (same "compile
            # third-party code as-is" rationale as -w -fno-strict-aliasing).
            if(NOT MSVC)
                set(_scotch_cflags -w -fno-strict-aliasing
                    -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0)
                target_compile_options(scotch    PRIVATE ${_scotch_cflags})
                target_compile_options(scotcherr PRIVATE ${_scotch_cflags})
                target_compile_options(esmumps   PRIVATE ${_scotch_cflags})
            endif()
            # The C caller sites are remapped bare→_mumps by force-including
            # scotch_rename_mumps.h. The Fortran caller (ana_orderings_
            # wrappers_m.F / [dz]ana_aux_par.F) reaches Scotch through the
            # scotchf.h interface and cannot force-include a C header, so
            # rename the 6 sequential Fortran entry points it actually calls
            # via CPP defines (the .F are preprocessed). CONTEXT*/DGRAPH*
            # bindings it also names are PT-Scotch/absent in this sequential
            # build and stay in guarded-out branches. -ffixed-line-length-none
            # (already set) absorbs the +6-char expansion in fixed form.
            set(_scotch_f_renames
                SCOTCHFGRAPHINIT=SCOTCHFGRAPHINIT_mumps
                SCOTCHFGRAPHBUILD=SCOTCHFGRAPHBUILD_mumps
                SCOTCHFGRAPHEXIT=SCOTCHFGRAPHEXIT_mumps
                SCOTCHFGRAPHPART=SCOTCHFGRAPHPART_mumps
                SCOTCHFSTRATINIT=SCOTCHFSTRATINIT_mumps
                SCOTCHFSTRATEXIT=SCOTCHFSTRATEXIT_mumps)
            set(MUMPS_HAVE_SCOTCH TRUE)
            message(STATUS "MUMPS: Scotch ordering enabled (ICNTL(7)=3)")
        endif()
    endif()
    if(NOT MUMPS_HAVE_SCOTCH)
        message(STATUS "MUMPS: Scotch ordering unavailable "
                       "(no _mumps_scotch_libsrc) — ICNTL(7)=3 will be inactive")
    endif()

    if(TARGET ${LIB_PREFIX}mumps)
        foreach(_dep ${LIB_PREFIX}scalapack ${LIB_PREFIX}lapack ${LIB_PREFIX}blas)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}mumps PUBLIC ${_dep})
            endif()
        endforeach()
        if(MPI_Fortran_FOUND)
            # PRIVATE so consumers can substitute libmpiseq at link time
            # without dragging libmpifort.so along the PUBLIC dep chain.
            target_link_libraries(${LIB_PREFIX}mumps PRIVATE MPI::MPI_Fortran)
        endif()
        # Multifloats: the migrator rewrites MPI_DOUBLE_PRECISION →
        # MPI_FLOAT64X2 (and reduction op tokens) in MUMPS Fortran sources,
        # then injects ``USE multifloats_mpi_f`` so those names resolve.
        # Link the module library so the .mod is on the include path and
        # the C++ side that actually defines the handles is pulled in.
        if(NEEDS_MULTIFLOATS AND TARGET multifloats_mpi_f)
            target_link_libraries(${LIB_PREFIX}mumps PUBLIC multifloats_mpi_f)
        endif()
        # Quad (kind16): same story — the migrator rewrites the reduce-op
        # tokens to MPI_QQ_SUM / MPI_XX_* on the standard MPI_REAL16 /
        # MPI_COMPLEX32 datatypes and injects ``USE quad_mpi_f``. Link the
        # shim so its .mod is on the include path and quad_mpi (which
        # defines the handles) is pulled in transitively.
        if(NEEDS_QUAD_MPI AND TARGET quad_mpi_f)
            target_link_libraries(${LIB_PREFIX}mumps PUBLIC quad_mpi_f)
        endif()
        # -Dpord reaches the migrated Fortran (CMAKE_Fortran_PREPROCESS is
        # ON) so [ey|qx|mw]ana_aux.F et al. compile the PORD call path
        # (MUMPS_PORDF, resolved from mumps_common's mumps_pord.c).
        if(MUMPS_HAVE_PORD)
            target_compile_definitions(${LIB_PREFIX}mumps PRIVATE pord)
        endif()
        # -Dmetis reaches the migrated analysis Fortran (ana_set_ordering.F,
        # ana_orderings_wrappers_m.F, [ey|qx|mw]ana_aux.F) so the METIS
        # nested-dissection call path (MUMPS_METIS_NODEND* → the namespaced
        # METIS_MUMPS_NODEND, resolved from mumps_common's mumps_metis.c)
        # compiles; without it those branches drop to the "ordering
        # unavailable" error.
        if(MUMPS_HAVE_METIS)
            target_compile_definitions(${LIB_PREFIX}mumps PRIVATE metis)
        endif()
        # -Dscotch reaches the migrated analysis Fortran (ana_set_ordering.F,
        # ana_orderings_wrappers_m.F, [ey|qx|mw]ana_aux.F) so the Scotch
        # ordering call path (MUMPS_SCOTCH* → the namespaced SCOTCH_*_mumps,
        # resolved from mumps_common's mumps_scotch.c) compiles; without
        # it those branches drop to the "ordering unavailable" error.
        if(MUMPS_HAVE_SCOTCH)
            target_compile_definitions(${LIB_PREFIX}mumps PRIVATE scotch
                $<$<COMPILE_LANGUAGE:Fortran>:${_scotch_f_renames}>)
            # [ey|qx|mw]ana_aux_par.F INCLUDE 'scotchf.h' (the Fortran
            # Scotch header) under -Dscotch — expose the Scotch include dir.
            target_include_directories(${LIB_PREFIX}mumps PRIVATE
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
        endif()
    endif()
    if(TARGET mumps_common AND MPI_Fortran_FOUND)
        target_link_libraries(mumps_common PRIVATE MPI::MPI_Fortran)
    endif()
    # ana_set_ordering.F / ana_orderings_wrappers_m.F / tools_common.F /
    # mumps_print_defined.F (in mumps_common) gate their PORD branches on
    # the same define.
    if(TARGET mumps_common AND MUMPS_HAVE_PORD)
        target_compile_definitions(mumps_common PRIVATE pord)
    endif()
    # Same for METIS: ana_set_ordering.F / ana_orderings_wrappers_m.F /
    # tools_common.F / mumps_print_defined.F gate their METIS branches on
    # the -Dmetis define too.
    if(TARGET mumps_common AND MUMPS_HAVE_METIS)
        target_compile_definitions(mumps_common PRIVATE metis)
    endif()
    # Same for Scotch: ana_set_ordering.F / ana_orderings_wrappers_m.F /
    # tools_common.F / mumps_print_defined.F gate their Scotch branches on
    # the -Dscotch define too.
    if(TARGET mumps_common AND MUMPS_HAVE_SCOTCH)
        target_compile_definitions(mumps_common PRIVATE scotch
            $<$<COMPILE_LANGUAGE:Fortran>:${_scotch_f_renames}>)
        # ana_orderings_wrappers_m.F INCLUDE 'scotchf.h' under -Dscotch.
        target_include_directories(mumps_common PRIVATE
            $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
    endif()

    # ── MUMPS C-side bridge ─────────────────────────────────────────
    # Builds ``${LIB_PREFIX}mumps_c`` (e.g. ``qmumps_c``), a packaged
    # C entry-point library that lets C consumers call MUMPS through
    # the migrated qmumps Fortran. Provides both real (``qmumps_c``)
    # and complex (``xmumps_c``) entry points by compiling upstream
    # ``mumps_c.c`` twice with per-arithmetic macro renames, plus
    # the upstream C runtime (mumps_common.c, IO, save/restore,
    # thread, numa, flytes, config, metis, pord, scotch).
    #
    # Per-target headers (``qmumps_c.h``, ``xmumps_c.h``,
    # ``mumps_c_types_extended.h``) come from
    # ``tests/mumps/target_${TARGET_NAME}/c/include/`` — produced by
    # the migrator at stage time. ``mumps_c_types_extended.h``
    # #includes upstream ``mumps_c_types.h`` and overrides its
    # ``double`` widths with the target's extended types; it is
    # force-included into the bridge objects below.
    #
    # The same machinery used to live inside tests/mumps/CMakeLists.txt
    # as the test-only ``mumps_c_bridge`` target; promoting it here
    # lets ``find_package(qmumps)`` deliver a usable C API instead of
    # the misleading raw ``dmumps_c.h`` (which expects libdmumps, not
    # libqmumps).
    set(_mumps_c_src ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_src)
    set(_mumps_c_inc ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_include)
    set(_mumps_c_bridge_inc_target
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/target_${TARGET_NAME}/c/include)
    set(_mumps_c_bridge_inc_shared
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/c/include)
    if(TARGET ${LIB_PREFIX}mumps
            AND IS_DIRECTORY ${_mumps_c_src}
            AND IS_DIRECTORY ${_mumps_c_inc}
            AND IS_DIRECTORY ${_mumps_c_bridge_inc_target})
        string(TOUPPER "${LIB_PREFIX}" _pfx_R)
        string(TOUPPER "${LIB_PREFIX_COMPLEX}" _pfx_C)
        # Per-arithmetic OBJECT libraries — compile ``mumps_c.c`` twice
        # with per-arith macro renames; the ``mumps_c.h`` ↔
        # ``dmumps_c.h``-style override gives us ``qmumps_c`` (real)
        # and ``xmumps_c`` (complex) entry points landing in the same
        # archive.
        add_library(${LIB_PREFIX}mumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
        target_include_directories(${LIB_PREFIX}mumps_c_obj PRIVATE
            ${_mumps_c_bridge_inc_target} ${_mumps_c_bridge_inc_shared}
            ${_mumps_c_inc} ${_mumps_c_src})
        target_compile_definitions(${LIB_PREFIX}mumps_c_obj PRIVATE
            MUMPS_ARITH=MUMPS_ARITH_d Add_
            dmumps_f77_=${LIB_PREFIX}mumps_f77_
            dmumps_set_tmp_ptr_=${LIB_PREFIX}mumps_set_tmp_ptr_
            dmumps_c=${LIB_PREFIX}mumps_c
            DMUMPS_STRUC_C=${_pfx_R}MUMPS_STRUC_C)
        target_compile_options(${LIB_PREFIX}mumps_c_obj PRIVATE
            -include mumps_c_types_extended.h)

        add_library(${LIB_PREFIX_COMPLEX}mumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
        target_include_directories(${LIB_PREFIX_COMPLEX}mumps_c_obj PRIVATE
            ${_mumps_c_bridge_inc_target} ${_mumps_c_bridge_inc_shared}
            ${_mumps_c_inc} ${_mumps_c_src})
        target_compile_definitions(${LIB_PREFIX_COMPLEX}mumps_c_obj PRIVATE
            MUMPS_ARITH=MUMPS_ARITH_z Add_
            zmumps_f77_=${LIB_PREFIX_COMPLEX}mumps_f77_
            zmumps_set_tmp_ptr_=${LIB_PREFIX_COMPLEX}mumps_set_tmp_ptr_
            zmumps_c=${LIB_PREFIX_COMPLEX}mumps_c
            ZMUMPS_STRUC_C=${_pfx_C}MUMPS_STRUC_C)
        target_compile_options(${LIB_PREFIX_COMPLEX}mumps_c_obj PRIVATE
            -include mumps_c_types_extended.h)

        # ── Shared, arithmetic-independent C runtime ────────────────
        # The upstream C runtime (mumps_common.c, IO, save/restore,
        # thread, numa, flytes, config, metis, pord, scotch) plus
        # mumps_addr.c carry ZERO arithmetic/width-macro references, so
        # a single compilation serves every arithmetic. They are folded
        # DIRECTLY INTO the shared Fortran common archive ``mumps_common``
        # — the C runtime's exact twin: also arithmetic-independent, also
        # the single shared archive that BOTH the migrated bridge
        # (${LIB_PREFIX}mumps_c) and the genuine double bridge
        # (dmumps_c / zmumps_c, below) already link. One merged archive
        # keeps exactly one copy of the runtime objects — the two
        # precisions still coexist in a single executable with no
        # duplicate-symbol clash — and drops the separate installed
        # ``mumps_c_runtime`` archive. The C sources compile byte-for-byte
        # as before (their include order + force-include are preserved
        # via C-language-gated per-target options below).
        #
        # mumps_addr.c stays at upstream double precision (the migrated
        # Fortran's xtools.F / qtools.F UPDATE_PROGRESS pass DP per
        # keep-kind manifest, and the genuine dmumps Fortran is double
        # throughout — compiling with the extended types would corrupt
        # the 8-byte caller buffer). It is a separate OBJECT so it can
        # opt out of the extended force-include the other runtime .c
        # tolerate as a verified no-op; it lists the upstream include
        # dir first so its mumps_c_types.h stays plain ``double``.
        add_library(mumps_c_runtime_addr_obj OBJECT ${_mumps_c_src}/mumps_addr.c)
        target_include_directories(mumps_c_runtime_addr_obj PRIVATE
            ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
        target_compile_definitions(mumps_c_runtime_addr_obj PRIVATE Add_)

        # Fold the C runtime .c (plus the DP-pinned mumps_addr object)
        # into the Fortran common archive. mumps_common is guaranteed here
        # (${LIB_PREFIX}mumps PUBLIC-links it via add_migrated_fortran_library).
        set(_mumps_c_runtime_srcs
            ${_mumps_c_src}/mumps_common.c
            ${_mumps_c_src}/mumps_io.c
            ${_mumps_c_src}/mumps_io_basic.c
            ${_mumps_c_src}/mumps_io_err.c
            ${_mumps_c_src}/mumps_io_thread.c
            ${_mumps_c_src}/mumps_save_restore_C.c
            ${_mumps_c_src}/mumps_register_thread.c
            ${_mumps_c_src}/mumps_thread.c
            ${_mumps_c_src}/mumps_thread_affinity.c
            ${_mumps_c_src}/mumps_numa.c
            ${_mumps_c_src}/mumps_flytes.c
            ${_mumps_c_src}/mumps_config_file_C.c
            ${_mumps_c_src}/mumps_metis.c
            ${_mumps_c_src}/mumps_metis64.c
            ${_mumps_c_src}/mumps_metis_int.c
            ${_mumps_c_src}/mumps_pord.c
            ${_mumps_c_src}/mumps_scotch.c
            ${_mumps_c_src}/mumps_scotch64.c
            ${_mumps_c_src}/mumps_scotch_int.c)
        target_sources(mumps_common PRIVATE
            $<TARGET_OBJECTS:mumps_c_runtime_addr_obj>
            ${_mumps_c_runtime_srcs})
        # Include dirs the runtime .c need. Harmless to the Fortran
        # sources (headers, not modules), so left ungated. PUBLIC so
        # C consumers of the installed common see the bridge headers.
        target_include_directories(mumps_common
            PUBLIC
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_target}>
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_shared}>
                $<BUILD_INTERFACE:${_mumps_c_inc}>
                $<BUILD_INTERFACE:${_mumps_c_src}>
                $<INSTALL_INTERFACE:include/mumps>)
        # ``Add_`` (C name-mangling) and the extended force-include are
        # C-only — gate them so the Fortran common compiles unchanged.
        target_compile_options(mumps_common PRIVATE
            $<$<COMPILE_LANGUAGE:C>:-include>
            $<$<COMPILE_LANGUAGE:C>:mumps_c_types_extended.h>)
        target_compile_definitions(mumps_common PRIVATE
            $<$<COMPILE_LANGUAGE:C>:Add_>)
        # PORD/METIS/SCOTCH ``-D`` defines are already on mumps_common
        # (set above for the Fortran ana_*/print_defined branches) and
        # apply to the C runtime too. Add the C-side include dirs (PORD
        # <space.h>, METIS <metis.h>) and PUBLIC-link the algorithm
        # archives so the ordering symbols mumps_pord.c / mumps_metis.c /
        # mumps_scotch.c reference resolve and are pulled through the
        # RESCAN umbrellas into the final link.
        if(MUMPS_HAVE_PORD)
            target_include_directories(mumps_common PRIVATE
                $<BUILD_INTERFACE:${_mumps_pord_inc}>)
            target_link_libraries(mumps_common PUBLIC pord)
        endif()
        if(MUMPS_HAVE_METIS)
            target_include_directories(mumps_common PRIVATE
                $<BUILD_INTERFACE:${_mumps_metis_inc}>)
            target_link_libraries(mumps_common PUBLIC metis)
        endif()
        # Scotch: mumps_scotch*.c reference *bare* SCOTCH_* names, so
        # force-include scotch_rename_mumps.h (C-only) to remap them onto
        # the suffixed archive symbols. The scotch include dir is already
        # on mumps_common (Fortran scotchf.h). SHELL keeps -include + path
        # as one group so CMake doesn't dedup it against the extended
        # force-include and drop the header to a stray input file.
        if(MUMPS_HAVE_SCOTCH)
            target_compile_options(mumps_common PRIVATE
                "SHELL:$<$<COMPILE_LANGUAGE:C>:-include ${_mumps_scotch_inc}/scotch_rename_mumps.h>")
            target_link_libraries(mumps_common PUBLIC esmumps scotch scotcherr)
        endif()
        if(MPI_C_FOUND)
            target_link_libraries(mumps_common PRIVATE MPI::MPI_C)
            target_link_libraries(mumps_common INTERFACE
                $<LINK_ONLY:MPI::MPI_C>)
        endif()

        # ── Fold the migrated C bridge INTO the Fortran archive ─────
        # The two per-arithmetic bridge objects are compiled INTO
        # ${LIB_PREFIX}mumps, so libeymumps.a carries both the e/y Fortran
        # solvers and their C entry points (${LIB_PREFIX}mumps_c() /
        # ${LIB_PREFIX_COMPLEX}mumps_c()) — no separate libeymumps_c.a.
        # The shared runtime now lives inside mumps_common (folded above),
        # which ${LIB_PREFIX}mumps already PUBLIC-links — so it is LINKED,
        # not bundled here, and isn't duplicated against the genuine bridge
        # in a combined link. The bridge↔solver reference (qmumps_c →
        # qmumps_f77_) is now internal to the one archive; only the
        # solver↔common (runtime) cycle remains, which the
        # ${LIB_PREFIX}mumps_full RESCAN below still resolves.
        target_sources(${LIB_PREFIX}mumps PRIVATE
            $<TARGET_OBJECTS:${LIB_PREFIX}mumps_c_obj>
            $<TARGET_OBJECTS:${LIB_PREFIX_COMPLEX}mumps_c_obj>)
        target_include_directories(${LIB_PREFIX}mumps
            PUBLIC
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_target}>
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_shared}>
                $<BUILD_INTERFACE:${_mumps_c_inc}>
                $<BUILD_INTERFACE:${_mumps_c_src}>
                $<INSTALL_INTERFACE:include/mumps>)
        if(MPI_C_FOUND)
            target_link_libraries(${LIB_PREFIX}mumps PRIVATE MPI::MPI_C)
            target_link_libraries(${LIB_PREFIX}mumps INTERFACE
                $<LINK_ONLY:MPI::MPI_C>)
        endif()

        # Umbrella INTERFACE target. The qmumps archive (Fortran solvers
        # + the folded-in C bridge) and mumps_common (Fortran common +
        # the folded-in C runtime) are mutually recursive at link time:
        # qmumps.a references mumps_size_c_ (in mumps_common's
        # mumps_addr.c) and the common is scanned back against qmumps.a.
        # A plain left-to-right scan misses one or the other. The
        # ``${LIB_PREFIX}mumps_full`` target wraps the archives in a
        # ``LINK_GROUP:RESCAN`` so ld --start-group / --end-group resolves
        # the cycle automatically. Consumers who use
        # ``kind16_libraries::${LIB_PREFIX}mumps_full`` don't have to know
        # about the static-archive ordering subtlety. Requires CMake 3.24+
        # on the consumer side. ``mumps_common`` is named explicitly in
        # the group because mumps_size_c_ (Fortran → runtime) now lives
        # there after the runtime was folded into the common archive.
        add_library(${LIB_PREFIX}mumps_full INTERFACE)
        # Propagate include directories (mumps_full is a meta-target;
        # by itself it has no headers, so pull them from qmumps
        # which now carries the per-target bridge include path).
        target_include_directories(${LIB_PREFIX}mumps_full INTERFACE
            $<INSTALL_INTERFACE:include/mumps>)
        # ``install(EXPORT)`` rewrites bare target names to
        # ``${PROJECT_NAME}::name`` in INTERFACE_LINK_LIBRARIES, but
        # only for top-level items — it doesn't recurse into the
        # comma-separated arg list of a ``LINK_GROUP`` genex. We
        # therefore embed the namespace literals here so the
        # serialized Targets file points at valid imported targets at
        # consumer time.
        target_link_libraries(${LIB_PREFIX}mumps_full INTERFACE
            "$<LINK_GROUP:RESCAN,${PROJECT_NAME}::${LIB_PREFIX}mumps,${PROJECT_NAME}::mumps_common>")

        # ══ Genuine double-precision MUMPS (additive) ═══════════════
        # Build real libdmumps / libzmumps + a plain-``double`` dmumps_c
        # / zmumps_c C bridge from the PRISTINE upstream sources staged
        # verbatim in _mumps_upstream_src, so the already-installed
        # upstream dmumps_c.h / zmumps_c.h headers finally have a
        # matching library behind them. Strictly additive: no migrated
        # target's emitted object code changes.
        #
        # Links system double-precision BLAS / LAPACK / ScaLAPACK
        # (find_package / find_library — NOT rebuilt here). If any is
        # missing the genuine build is skipped with a STATUS message so
        # the extended-precision release still configures cleanly.
        #
        # NOTE: genuine s/c/d/z and the extended-precision solvers now
        # share ONE common — the migrated ``mumps_common``. Its 36 files
        # are pinned verbatim from upstream (recipe ``copy_files``), so
        # they carry no working precision and serve every arithmetic
        # identically — exactly as upstream compiles ``libmumps_common``
        # once and links it into libsmumps/libdmumps/libcmumps/libzmumps
        # alike. This retires the former separate ``dzmumps_common``: with
        # every common module writing into the single shared Fortran
        # module directory, one producer keeps the .mod set unambiguous
        # for both the genuine and the extended solvers.
        find_package(BLAS QUIET)
        find_package(LAPACK QUIET)
        # ── Parallel stack for the genuine solvers ──────────────────
        # Prefer the in-tree double-precision parallel stack (ScaLAPACK
        # / PBLAS / BLACS / LAPACK / BLAS built above from the SAME MPI
        # this release configured with). Being compiled against the
        # chosen MPI it is ABI-compatible with that MPI by construction,
        # so the genuine build no longer hinges on a SYSTEM ScaLAPACK
        # whose BLACS flavor happens to match the active MPI — the part
        # that is fragile under, e.g., Intel MPI. Each std archive (plain
        # name, from add_standard_*) is paired with its migrated
        # ``_common`` companion (from add_migrated_*); keep whichever
        # candidates are real targets in this configuration. This is the
        # same archive set the scratchpad/iso genuine-``gen`` link uses.
        set(_bundled_parallel_candidates
            scalapack scalapack_common
            pblas pblas_common pbblas
            ptzblas ptzblas_common
            blacs blacs_common
            lapack lapack_common blas)
        set(_genuine_parallel_stack "")
        foreach(_t ${_bundled_parallel_candidates})
            if(TARGET ${_t})
                list(APPEND _genuine_parallel_stack ${_t})
            endif()
        endforeach()
        if(TARGET scalapack AND TARGET pblas AND TARGET blacs
                AND TARGET lapack AND TARGET blas)
            set(_have_bundled_parallel TRUE)
        else()
            set(_have_bundled_parallel FALSE)
        endif()

        # System double-precision ScaLAPACK — the fallback when the
        # bundled stack is absent (e.g. BASELINE_BUILD emits no migrated
        # companions and skips the std parallel archives). ScaLAPACK
        # ships no CMake package and is MPI-flavored on Debian
        # (libscalapack-openmpi / libscalapack-mpich): resolve the
        # variant matching the MPI flavor already selected above; allow
        # an override via -DMUMPS_SCALAPACK_LIBRARY=<path>.
        if(NOT DEFINED MUMPS_SCALAPACK_LIBRARY OR NOT MUMPS_SCALAPACK_LIBRARY)
            set(_scalapack_names scalapack)
            if(MPI_TAG MATCHES "openmpi")
                list(INSERT _scalapack_names 0 scalapack-openmpi)
            elseif(MPI_TAG MATCHES "mpich")
                list(INSERT _scalapack_names 0 scalapack-mpich)
            endif()
            find_library(MUMPS_SCALAPACK_LIBRARY NAMES ${_scalapack_names})
        endif()
        set(_scalapack_lib "${MUMPS_SCALAPACK_LIBRARY}")
        if(_have_bundled_parallel)
            set(_have_system_parallel FALSE)  # bundled stack wins
        elseif(LAPACK_FOUND AND BLAS_FOUND AND _scalapack_lib)
            set(_have_system_parallel TRUE)
        else()
            set(_have_system_parallel FALSE)
        endif()

        if((_have_bundled_parallel OR _have_system_parallel)
                AND MPI_Fortran_FOUND
                AND IS_DIRECTORY ${_mumps_c_src}
                AND IS_DIRECTORY ${_mumps_c_inc})
            # ── source lists ─────────────────────────────────────────
            # d-arith = glob(d*.F) minus the sole d-prefixed common file
            # double_linked_list.F; z-arith = glob(z*.F) (no z-common).
            file(GLOB _dmumps_F CONFIGURE_DEPENDS ${_mumps_c_src}/d*.F)
            list(FILTER _dmumps_F EXCLUDE REGEX "/double_linked_list\\.F$")
            file(GLOB _zmumps_F CONFIGURE_DEPENDS ${_mumps_c_src}/z*.F)
            # s-arith = glob(s*.F) minus the three unprefixed ``sol_*`` commons
            # that also match s* and live in the shared ``mumps_common`` (the
            # ssol_*/sfac_sol_*/smumps_sol_* files are genuine single solver
            # sources and stay). c-arith = glob(c*.F) (no c-prefixed common).
            file(GLOB _smumps_F CONFIGURE_DEPENDS ${_mumps_c_src}/s*.F)
            list(FILTER _smumps_F EXCLUDE REGEX "/sol_(common|ds_common_m|omp_common_m)\\.F$")
            file(GLOB _cmumps_F CONFIGURE_DEPENDS ${_mumps_c_src}/c*.F)
            # The arithmetic-independent common is NOT rebuilt here: the
            # migrated ``mumps_common`` already holds all 36 upstream common
            # files verbatim (recipe ``copy_files``) and is the sole producer
            # of their .mod set. The genuine solvers link it below, exactly
            # as the extended ${LIB_PREFIX}mumps does.

            # ── genuine solvers, combined by real+complex family ─────
            # Mirrors the extended one-archive layout (libeymumps.a holds
            # both the e-real and y-complex solvers): the genuine double
            # real (d) and complex (z) sources compile into ONE archive
            # libdzmumps.a, and the single real (s) + complex (c) into
            # libscmumps.a. Upstream module namespaces are arith-prefixed
            # (DMUMPS_*/ZMUMPS_*, SMUMPS_*/CMUMPS_*) with zero overlap, so
            # a family's two arithmetics coexist in one target cleanly.
            # ${ARITH}mumps_gpu.c is the sole C per arithmetic (the Fortran
            # calls its ${s,d,c,z}mumps_gpu_return symbol).
            add_library(dzmumps STATIC
                ${_dmumps_F} ${_mumps_c_src}/dmumps_gpu.c
                ${_zmumps_F} ${_mumps_c_src}/zmumps_gpu.c)
            fortran_module_layout(dzmumps)
            add_library(scmumps STATIC
                ${_smumps_F} ${_mumps_c_src}/smumps_gpu.c
                ${_cmumps_F} ${_mumps_c_src}/cmumps_gpu.c)
            fortran_module_layout(scmumps)

            foreach(_g dzmumps scmumps)
                # _mumps_c_bridge_inc_shared carries mumps_int_def.h,
                # which the gpu-stub C pulls in via mumps_c_types.h.
                target_include_directories(${_g} PRIVATE
                    ${_mumps_c_src} ${_mumps_c_inc} ${_mumps_c_bridge_inc_shared})
                # PRIVATE MPI Fortran on every genuine Fortran archive —
                # the sources #include mpif.h, and PRIVATE keeps the
                # libmpiseq substitution open (mirrors the migrated
                # mumps_common / ${LIB_PREFIX}mumps).
                target_link_libraries(${_g} PRIVATE MPI::MPI_Fortran)
                # Add_ (C name-mangling) only for the gpu-stub C; don't
                # pass -DAdd_ to the Fortran compile.
                target_compile_definitions(${_g} PRIVATE
                    $<$<COMPILE_LANGUAGE:C>:Add_>)
                # -Dpord for the genuine analysis Fortran ([dz]ana_aux.F,
                # ana_set_ordering.F, …); MUMPS_PORDF resolves from the
                # shared mumps_common (folded C runtime), which links libpord.
                if(MUMPS_HAVE_PORD)
                    target_compile_definitions(${_g} PRIVATE
                        $<$<COMPILE_LANGUAGE:Fortran>:pord>)
                endif()
                # -Dmetis for the genuine analysis Fortran ([dz]ana_aux.F,
                # ana_set_ordering.F, …); the METIS_MUMPS_* entry points
                # resolve from the shared mumps_common (folded C runtime),
                # which links libmetis.
                if(MUMPS_HAVE_METIS)
                    target_compile_definitions(${_g} PRIVATE
                        $<$<COMPILE_LANGUAGE:Fortran>:metis>)
                endif()
                # -Dscotch for the genuine analysis Fortran ([dz]ana_aux.F,
                # ana_set_ordering.F, …); the SCOTCH_*_mumps entry points
                # resolve from the shared mumps_common (folded C runtime),
                # which links scotch/esmumps.
                if(MUMPS_HAVE_SCOTCH)
                    target_compile_definitions(${_g} PRIVATE
                        $<$<COMPILE_LANGUAGE:Fortran>:scotch>
                        $<$<COMPILE_LANGUAGE:Fortran>:${_scotch_f_renames}>)
                    # [dz]ana_aux_par.F / ana_orderings_wrappers_m.F
                    # INCLUDE 'scotchf.h' under -Dscotch.
                    target_include_directories(${_g} PRIVATE
                        $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
                endif()
                # Pristine upstream compiles within column limits, but
                # relieve line length anyway to match the migrated macro
                # and stay robust to compiler defaults.
                if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran"
                   OR FORTRAN_COMPILER_FAMILY STREQUAL "flang")
                    target_compile_options(${_g} PRIVATE
                        $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>
                        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
                elseif(FORTRAN_COMPILER_FAMILY STREQUAL "intel")
                    target_compile_options(${_g} PRIVATE
                        $<$<COMPILE_LANGUAGE:Fortran>:-extend-source>)
                endif()
            endforeach()

            # Solvers link the genuine common (propagates its .mod dir),
            # then the parallel stack + MPI Fortran. Every parallel dep
            # is wrapped in $<BUILD_INTERFACE:> so it builds our tests but
            # is NOT forced into install(EXPORT) — downstream consumers
            # list their own BLAS/LAPACK/ScaLAPACK on the link line,
            # matching standard MUMPS practice and the repo convention
            # that no cross-package find_dependency is emitted. (The
            # dzmumps_full/scmumps_full umbrellas below RESCAN only the
            # solver + bridge + common (which now carries the C runtime),
            # never the parallel stack, so this stays a pure build-time
            # dependency.)
            foreach(_s dzmumps scmumps)
                target_link_libraries(${_s} PUBLIC mumps_common)
                if(_have_bundled_parallel)
                    # MPI-matched in-tree stack. The circular refs between
                    # these static archives are resolved by the RESCAN
                    # link group each test builds (tests/mumps); listing
                    # them PUBLIC here mirrors the migrated ${LIB_PREFIX}mumps.
                    foreach(_p ${_genuine_parallel_stack})
                        target_link_libraries(${_s} PUBLIC $<BUILD_INTERFACE:${_p}>)
                    endforeach()
                else()
                    target_link_libraries(${_s} PUBLIC
                        $<BUILD_INTERFACE:LAPACK::LAPACK>
                        $<BUILD_INTERFACE:BLAS::BLAS>
                        $<BUILD_INTERFACE:${_scalapack_lib}>)
                endif()
            endforeach()

            # ── genuine C bridges (native precision, real names) ─────
            # mumps_c.c compiled per-arith with NO rename and NO extended
            # header → genuine s/d/c/z mumps_c() at upstream native
            # precision, backing the installed ?mumps_c.h headers.
            # Upstream include dir FIRST so mumps_c_types.h stays at the
            # native float/double widths.
            add_library(dmumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
            target_include_directories(dmumps_c_obj PRIVATE
                ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
            target_compile_definitions(dmumps_c_obj PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_d Add_)
            add_library(zmumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
            target_include_directories(zmumps_c_obj PRIVATE
                ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
            target_compile_definitions(zmumps_c_obj PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_z Add_)
            add_library(smumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
            target_include_directories(smumps_c_obj PRIVATE
                ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
            target_compile_definitions(smumps_c_obj PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_s Add_)
            add_library(cmumps_c_obj OBJECT ${_mumps_c_src}/mumps_c.c)
            target_include_directories(cmumps_c_obj PRIVATE
                ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
            target_compile_definitions(cmumps_c_obj PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_c Add_)

            # Fold the per-family C bridge objects INTO the solver
            # archives: libdzmumps.a carries both the d/z Fortran solvers
            # and their C entry points (dmumps_c()+zmumps_c()), libscmumps.a
            # the s/c pair — no separate libdzmumps_c.a / libscmumps_c.a,
            # matching the folded extended ${LIB_PREFIX}mumps. The shared
            # runtime now lives inside mumps_common (linked below), not
            # bundled. The bridge↔solver reference is now internal to the
            # archive; the solver↔common (runtime) cycle is resolved by the
            # *_full RESCAN below.
            target_sources(dzmumps PRIVATE
                $<TARGET_OBJECTS:dmumps_c_obj> $<TARGET_OBJECTS:zmumps_c_obj>)
            target_sources(scmumps PRIVATE
                $<TARGET_OBJECTS:smumps_c_obj> $<TARGET_OBJECTS:cmumps_c_obj>)
            foreach(_g dzmumps scmumps)
                target_include_directories(${_g}
                    PUBLIC
                        $<BUILD_INTERFACE:${_mumps_c_inc}>
                        $<BUILD_INTERFACE:${_mumps_c_src}>
                        $<INSTALL_INTERFACE:include/mumps>)
                # The shared C runtime is reached via mumps_common (linked
                # PUBLIC in the solver loop above); nothing extra here.
                if(MPI_C_FOUND)
                    target_link_libraries(${_g} PRIVATE MPI::MPI_C)
                    target_link_libraries(${_g} INTERFACE
                        $<LINK_ONLY:MPI::MPI_C>)
                endif()
            endforeach()

            # Umbrella INTERFACE targets — same solver↔common archive cycle
            # as ${LIB_PREFIX}mumps_full; RESCAN over the (combined) solver
            # with its folded-in bridge and mumps_common (Fortran common +
            # the folded-in shared C runtime).
            add_library(dzmumps_full INTERFACE)
            target_include_directories(dzmumps_full INTERFACE
                $<INSTALL_INTERFACE:include/mumps>)
            target_link_libraries(dzmumps_full INTERFACE
                "$<LINK_GROUP:RESCAN,${PROJECT_NAME}::dzmumps,${PROJECT_NAME}::mumps_common>")
            add_library(scmumps_full INTERFACE)
            target_include_directories(scmumps_full INTERFACE
                $<INSTALL_INTERFACE:include/mumps>)
            target_link_libraries(scmumps_full INTERFACE
                "$<LINK_GROUP:RESCAN,${PROJECT_NAME}::scmumps,${PROJECT_NAME}::mumps_common>")

            if(_have_bundled_parallel)
                message(STATUS
                    "MUMPS: genuine single+double enabled (libdzmumps + libscmumps; MPI-matched in-tree ScaLAPACK stack)")
            else()
                message(STATUS
                    "MUMPS: genuine single+double enabled (libdzmumps + libscmumps; system ScaLAPACK=${_scalapack_lib})")
            endif()
        else()
            message(STATUS
                "MUMPS: genuine single+double build skipped — need the in-tree parallel stack OR system BLAS/LAPACK/ScaLAPACK, plus MPI Fortran "
                "(bundled=${_have_bundled_parallel} BLAS_FOUND=${BLAS_FOUND} LAPACK_FOUND=${LAPACK_FOUND} system_ScaLAPACK=${_scalapack_lib} MPI_Fortran=${MPI_Fortran_FOUND})")
        endif()
    endif()
endif()

# 10. libmpiseq — sequential MPI stub from MUMPS 5.8.2 upstream.
#
# Provides single-process implementations of the MPI / BLACS /
# ScaLAPACK symbols the migrated qmumps archive calls at link time
# (``mpi_init_``, ``mpi_send_``, ``mpi_bcast_``, ``blacs_pinfo_``,
# ``descinit_``, ``ddot_``, ...). Linking a test executable against
# ``mpiseq`` instead of ``MPI::MPI_Fortran`` (and against a stock
# BLACS/PBLAS/ScaLAPACK link line) lets it run as a plain executable,
# no ``mpiexec`` needed, even on machines without MPI installed.
#
# Constraints: NPROCS=1 only — collectives that require multi-rank
# coordination print "should not be called" and abort.
#
# ABI matching with the configured MPI: libseq ships its own mpif.h /
# mpi.h with hardcoded constant values (MPI_DOUBLE_PRECISION = 12,
# MPI_INTEGER = 13, ...) and its own COMMON-block names
# (/mpif_libseq/MPI_IN_PLACE). Migrated archives are compiled against
# the configured MPI's headers (Intel MPI: MPI_DOUBLE_PRECISION =
# 1275070495, MPI_IN_PLACE in /MPIPRIV1/, ...). For libseq's stubs to
# correctly identify the in-place sentinel and dispatch on datatypes
# in MUMPS_COPY, the stubs must see the SAME constants as the callers.
# We achieve this by overlaying the configured MPI's mpif.h / mpi.h
# onto _mpiseq_src/ at configure time — libseq's source still uses
# ``INCLUDE 'mpif.h'`` / ``#include "mpi.h"`` literally, but those
# resolve to the real-MPI copies. Skipped silently if find_package(MPI)
# didn't locate the headers.
#
# Sources are staged into ``_mpiseq_src/`` by ``migrator stage`` (see
# ``_std_dirs`` in ``src/migrator/__main__.py``); the build silently
# skips when that directory is absent (e.g., a bare repo with no
# upstream MUMPS vendored).
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src AND MPISEQ_IS_INTEL_MPI)
    # Replace libseq's mpif.h / mpi.h with the configured MPI's so
    # libmpiseq.a is ABI-compatible with archives compiled against
    # real MPI: same constant values (MPI_INTEGER = 1275069467, etc.)
    # and same COMMON-block layout (MPI_IN_PLACE in /MPIPRIV1/).
    # Implementation: delete libseq's bundled headers from the staging
    # tree and prepend the configured MPI's include dir to mpiseq's
    # PRIVATE include path. Fortran INCLUDE 'mpif.h' / C #include
    # "mpi.h" then resolve to the real-MPI copies, which also pulls in
    # auxiliary headers like mpio.h / mpiof.h that live in the same
    # directory. Skipped silently if find_package(MPI) didn't locate
    # the headers (libseq's bundled copies are then used as-is).
    if(MPI_Fortran_FOUND)
        file(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/mpif.h)
    endif()
    if(MPI_C_FOUND)
        file(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/mpi.h)
    endif()
    add_library(mpiseq STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/mpi.f
        # libseq's bundled mpic.c is replaced by cmake/mpiseq_c_stubs.c
        # (added below). Intel MPI defines MPI_Comm_f2c as a function-
        # like macro; libseq's mpic.c declares it as a function and the
        # two collide before any compile-flag tricks. Replacing the
        # whole file is cleaner.
        ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src/elapse.c
        # Per-precision ScaLAPACK stubs covering the Q/X/E/Y prefixes
        # the migrated qmumps/xmumps/emumps/ymumps archives request
        # (upstream's mpi.f only ships PD* / DOUBLE PRECISION variants).
        # Authored in-tree because external/ is read-only. M/W
        # (multifloats) live in mpiseq_mw_stubs.f90, appended below
        # only when NEEDS_MULTIFLOATS — they need TYPE(real64x2) /
        # TYPE(cmplx64x2) from the multifloats Fortran module.
        ${CMAKE_CURRENT_SOURCE_DIR}/mpiseq_qx_stubs.f
    )
    if(NEEDS_MULTIFLOATS AND TARGET multifloatsf)
        target_sources(mpiseq PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/mpiseq_mw_stubs.f90)
        target_link_libraries(mpiseq PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()
    # Supplementary C-side stubs (MPI_Send / MPI_Type_* / MPI_Comm_split
    # / MPI_Op_create / etc.) covering the C-MPI surface that BLACS /
    # PBLAS / multifloats_mpi reach for. Always folded in so a
    # downstream consumer can link this libmpiseq.a against the
    # archives in place of real MPI without further plumbing.
    target_sources(mpiseq PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/mpiseq_c_stubs.c)
    target_include_directories(mpiseq PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src>
    )
    # Prepend configured MPI's include path so libseq's INCLUDE 'mpif.h'
    # and #include "mpi.h" resolve to the real-MPI copies (we deleted
    # libseq's bundled headers above when MPI was found).
    if(MPI_Fortran_FOUND AND MPI_Fortran_INCLUDE_DIRS)
        target_include_directories(mpiseq BEFORE PRIVATE
            ${MPI_Fortran_INCLUDE_DIRS})
    endif()
    if(MPI_C_FOUND AND MPI_C_INCLUDE_DIRS)
        target_include_directories(mpiseq BEFORE PRIVATE
            ${MPI_C_INCLUDE_DIRS})
    endif()
    # Add_ matches gfortran's lowercase + trailing-underscore Fortran
    # symbol convention: elapse.h's #define mumps_elapse mumps_elapse_
    # fires under Add_, so the Fortran-side mpi.f's MPI_WTIME (which
    # calls ``mumps_elapse_``) resolves cleanly.
    target_compile_definitions(mpiseq PRIVATE Add_)
    # libseq's bundled mpic.c uses ``LIBSEQ_CALL`` (calling-convention
    # tag) and ``LIBSEQ_INT`` / ``LIBSEQ_INT8`` (Fortran-INTEGER shadow
    # types) macros that came from libseq's own mpi.h. With the
    # configured MPI's mpi.h overlaid those macros aren't defined —
    # plumb them in here. ``LIBSEQ_INT`` matches the Fortran INTEGER
    # default kind on this build (always ``int`` for our gfortran
    # target since we don't pass -fdefault-integer-8); ``LIBSEQ_INT8``
    # always picks the 8-byte kind for MUMPS_CHECKADDREQUAL's address
    # arithmetic.
    target_compile_definitions(mpiseq PRIVATE
        LIBSEQ_CALL=
        LIBSEQ_INT=int
        LIBSEQ_INT8=int64_t)
    fortran_module_layout(mpiseq)
endif()

# ── Install ──────────────────────────────────────────────────────────
# Each precision library gets its own export set and Config.cmake,
# so consumers can do: find_package(qblas), find_package(qlapack), etc.
# The common library is bundled into the same export as its precision library.

# Libraries that link MPI (transitively or directly) — their installs
# pick up an extra ``-${MPI_TAG}`` suffix so OpenMPI/MPICH/Intel-MPI
# variants can coexist in the same install prefix.
set(_MPI_DEPENDENT_LIBS blacs pblas pbblas scalapack mumps)


# install_library_headers(<lib_name>): install any ``*.h`` headers
# under ``${lib_name}/src/`` to ``include/${lib_name}/`` and extend
# every related target's INTERFACE include directories so consumers
# linking against the installed library pick up the include path via
# ``find_package``. Per-library subdirs are kept as a defensive
# convention against header-name clashes across libraries. The one
# historical clash — PBLAS and ScaLAPACK_C both shipping a ``pblas.h``
# with different contents — is now moot: ScaLAPACK_C's C sources fold
# into the ScaLAPACK archive and its internal headers (``pblas.h``,
# ``redist.h``, …) are never installed. No two installed libraries
# currently share a header basename.
#
# BLACS is excluded: ``Bdef.h`` and ``Bconfig.h`` are internal
# config headers (no user-facing C API ships with upstream BLACS —
# callers write their own ``extern void Cblacs_*`` declarations).
# Bdef.h additionally defines ``QCOMPLEX`` with field names that
# collide with xblas/blas_enum.h's ``QCOMPLEX`` when both are
# included in the same translation unit. Skipping the install avoids
# the clash without changing any source.
set(_HEADER_INSTALL_SKIP blacs)
function(install_library_headers lib_name)
    if("${lib_name}" IN_LIST _HEADER_INSTALL_SKIP)
        return()
    endif()
    set(_src ${CMAKE_CURRENT_SOURCE_DIR}/${lib_name}/src)
    file(GLOB _hdrs ${_src}/*.h)
    if(NOT _hdrs)
        return()
    endif()
    foreach(_t ${LIB_PREFIX}${lib_name} ${lib_name} ${lib_name}_common)
        if(TARGET ${_t})
            target_include_directories(${_t}
                INTERFACE $<INSTALL_INTERFACE:include/${lib_name}>)
        endif()
    endforeach()
    install(DIRECTORY ${_src}/
        DESTINATION include/${lib_name}
        FILES_MATCHING PATTERN "*.h")
endfunction()

# Compute the package name a standard archive is installed under.
# `blas` → `eplinalgStdBlas`, `scalapack` → `eplinalgStdScalapack`, etc.
# The capitalization keeps the bare lib name visually intact while
# satisfying CMake's package-name conventions (leading capital aids
# readability when grep'ing for find_package(eplinalg...) calls).
function(_std_pkg_name lib_name out_var)
    string(SUBSTRING "${lib_name}" 0 1 _first)
    string(SUBSTRING "${lib_name}" 1 -1 _rest)
    string(TOUPPER "${_first}" _first)
    set(${out_var} "eplinalgStd${_first}${_rest}" PARENT_SCOPE)
endfunction()

# Install the standard-precision archive ``lib_name`` (e.g. ``blas``,
# ``lapack``) into its own package ``eplinalgStd<LibName>``. The archive
# is byte-identical across per-precision build invocations (vanilla
# Netlib sources, no real-promotion flag, no precision-specific
# defines), so re-installing on top of an existing prefix overwrites
# matching content — which is exactly what consumers want: one
# ``eplinalg::blas`` target shared by ``find_package(qblas)``,
# ``find_package(eblas)``, ``find_package(mblas)``, ... avoiding the
# multi-target redefinition that bundling the std archive in each
# per-precision export used to cause.
function(_install_standard_archive lib_name)
    if(NOT TARGET ${lib_name})
        return()
    endif()
    _std_pkg_name(${lib_name} _pkg)
    set(_export "${_pkg}Targets")

    # In a given configure, multiple precision libs might trigger
    # install of the same std archive — guard so only the first wins.
    get_property(_done GLOBAL PROPERTY _FC_STD_INSTALLED_${_pkg})
    if(_done)
        return()
    endif()
    set_property(GLOBAL PROPERTY _FC_STD_INSTALLED_${_pkg} TRUE)

    set(_mpi_arg "")
    if("${lib_name}" IN_LIST _MPI_DEPENDENT_LIBS)
        set(_mpi_arg "MPI")
    endif()

    # Walk the std archive's INTERFACE_LINK_LIBRARIES for std-sibling
    # references (e.g. ``lapack`` PUBLIC-links ``blas``) so the
    # generated eplinalgStd<Lib>Config does find_dependency on the other
    # std packages it transitively needs.
    set(_std_deps "")
    get_target_property(_ifl ${lib_name} INTERFACE_LINK_LIBRARIES)
    if(_ifl)
        foreach(_d IN LISTS _ifl)
            if(TARGET "${_d}")
                _std_pkg_name("${_d}" _dep_pkg)
                if(NOT "${_dep_pkg}" STREQUAL "${_pkg}")
                    list(APPEND _std_deps "${_dep_pkg}")
                endif()
            endif()
        endforeach()
    endif()

    fortran_install_library(${lib_name}
        EXPORT ${_export}
        NAMESPACE eplinalg::
        ${_mpi_arg}
        DEPENDS ${_std_deps})
    get_target_property(_fc_moddir ${lib_name} Fortran_MODULE_DIRECTORY)
    if(_fc_moddir)
        fortran_install_modules(${lib_name})
    endif()
endfunction()

function(_install_library_pair lib_name)
    set(_precision_target "${LIB_PREFIX}${lib_name}")
    set(_common_target "${lib_name}_common")
    set(_export "${_precision_target}Targets")

    if(NOT TARGET ${_precision_target})
        return()
    endif()

    set(_install_mpi_arg "")
    if("${lib_name}" IN_LIST _MPI_DEPENDENT_LIBS)
        set(_install_mpi_arg "MPI")
        # MPI-dependent libs only compile when MPI is available; the
        # target exists but the build never produces an archive file,
        # which would later cause cmake --install to fail with "file
        # INSTALL cannot find ...". Skip install when MPI is missing —
        # the same gate that tests/blacs/pblas/scalapack subdirs apply.
        if(NOT MPI_C_FOUND OR NOT MPI_Fortran_FOUND)
            return()
        endif()
    endif()

    # Install per-library headers and extend the targets' install-time
    # include directories before the export captures them.
    install_library_headers(${lib_name})

    # Install the std archive into its own export FIRST. Required so
    # install(EXPORT ${_export}) below can resolve the precision
    # target's PUBLIC link to ``${lib_name}`` (the std archive) as a
    # cross-export reference ``eplinalg::${lib_name}`` — without this,
    # CMake fails with "target X not in any export set".
    _install_standard_archive(${lib_name})

    # If the std archive exists, its package is a transparent dep of
    # the precision Config so consumers only need find_package(qblas)
    # and ``eplinalg::blas`` resolves automatically.
    set(_precision_deps "")
    if(TARGET ${lib_name} AND NOT "${lib_name}" STREQUAL "${_precision_target}")
        _std_pkg_name(${lib_name} _std_pkg)
        list(APPEND _precision_deps "${_std_pkg}")
    endif()

    # NEEDS_MULTIFLOATS precision targets also PUBLIC-link the
    # la_constants_mw / la_xisnan_mw helper archives (see ~line 403),
    # each installed as its own Config package (EXPORT ${_t}Targets,
    # ~line 1439). Like eplinalgStdBlas they are transparent deps, so
    # consumers only need find_package(mblas) — without these the
    # generated Config references undefined imported targets and sets
    # mblas_FOUND = FALSE.
    if(NEEDS_MULTIFLOATS)
        foreach(_mf_helper la_constants_mw la_xisnan_mw)
            if(TARGET ${_mf_helper})
                list(APPEND _precision_deps ${_mf_helper})
            endif()
        endforeach()
    endif()

    # Install precision-specific library FIRST. fortran_install_library
    # only writes Config.cmake on the first call per export, scanning
    # THAT target's INTERFACE_LINK_LIBRARIES for cross-package deps.
    # The precision target (qmumps) has the canonical dependency set
    # for find_package consumers; siblings sharing the export
    # (qmumps_c, mumps_common) follow after.
    fortran_install_library(${_precision_target}
        EXPORT ${_export}
        NAMESPACE eplinalg::
        ${_install_mpi_arg}
        DEPENDS ${_precision_deps})
    get_target_property(_fc_moddir ${_precision_target} Fortran_MODULE_DIRECTORY)
    if(_fc_moddir)
        fortran_install_modules(${_precision_target})
    endif()

    # MUMPS-special: the C-side bridge is now folded into
    # ${LIB_PREFIX}mumps and the shared C runtime into mumps_common (both
    # installed via the generic paths), so the export just needs the
    # ordering leaves mumps_common PUBLIC-links and the umbrella INTERFACE
    # (${LIB_PREFIX}mumps_full) that wraps solver + common in a
    # LINK_GROUP:RESCAN. ``find_package(qmumps)`` then brings the Fortran
    # archive, its folded C entry points (qmumps_c, xmumps_c) and the
    # umbrella.
    if("${lib_name}" STREQUAL "mumps" AND TARGET mumps_common)
        # PORD ordering archive — PUBLIC-linked by mumps_common (folded C
        # runtime), so it must join the same export (CMake rejects an
        # exported target that requires a non-exported one). A leaf C
        # archive: no Fortran modules, no MPI, its include dir is
        # BUILD_INTERFACE-only, so a plain install(TARGETS) suffices.
        # Consumers of find_package(qmumps) need it on the link line to
        # resolve PORD's SPACE_ordering et al.
        if(TARGET pord)
            install(TARGETS pord EXPORT ${_export})
        endif()
        # METIS ordering archive — PUBLIC-linked by mumps_common, so it
        # must join the same export. Like pord, a leaf C archive with a
        # BUILD_INTERFACE-only include dir, so a plain install(TARGETS)
        # suffices; consumers need it to resolve METIS_MUMPS_NodeND et al.
        if(TARGET metis)
            install(TARGETS metis EXPORT ${_export})
        endif()
        # Scotch ordering archives — PUBLIC-linked by mumps_common
        # (esmumps + scotch + scotcherr), so they must join the same
        # export. Leaf C archives with BUILD_INTERFACE-only include dirs;
        # consumers need them to resolve SCOTCH_*_mumps / esmumps.
        if(TARGET scotch)
            install(TARGETS scotch scotcherr esmumps EXPORT ${_export})
        endif()
        install(TARGETS ${LIB_PREFIX}mumps_full
            EXPORT ${_export})

        # Genuine single+double solvers (additive; present only when the
        # system BLAS/LAPACK/ScaLAPACK deps were found). Each now carries
        # its folded-in C bridge, so there are no separate *_c archives.
        # dzmumps_common is retired — the genuine solvers link the shared
        # migrated ``mumps_common`` (installed via the library-pair path
        # below), so it is not listed here.
        foreach(_g dzmumps scmumps)
            if(TARGET ${_g})
                fortran_install_library(${_g} EXPORT ${_export} NAMESPACE eplinalg:: MPI)
                get_target_property(_g_moddir ${_g} Fortran_MODULE_DIRECTORY)
                if(_g_moddir)
                    fortran_install_modules(${_g})
                endif()
            endif()
        endforeach()
        foreach(_gf dzmumps_full scmumps_full)
            if(TARGET ${_gf})
                install(TARGETS ${_gf} EXPORT ${_export})
            endif()
        endforeach()
    endif()

    # Bundle common library into the same export
    if(TARGET ${_common_target})
        fortran_install_library(${_common_target} EXPORT ${_export} NAMESPACE eplinalg:: ${_install_mpi_arg})
        get_target_property(_fc_moddir ${_common_target} Fortran_MODULE_DIRECTORY)
        if(_fc_moddir)
            fortran_install_modules(${_common_target})
        endif()
    endif()
endfunction()

if(NOT BASELINE_BUILD)
    foreach(_lib ${STAGED_LIBRARIES})
        _install_library_pair(${_lib})
    endforeach()
endif()

# MUMPS additionally ships C-API entry-point headers in three layers:
#
#   1. ``_mumps_upstream_include/*.h`` — upstream's dmumps_c.h,
#      zmumps_c.h, dmumps_struc.h, mumps_compat.h and the plain-
#      ``double`` mumps_c_types.h (extended, not overwritten, below).
#   2. ``tests/mumps/c/include/mumps_int_def.h`` — fixed 32-bit
#      MUMPS_INT stub (upstream's is generated from a .in template).
#   3. ``tests/mumps/target_${TARGET_NAME}/c/include/*.h`` — the
#      standalone wrappers ``qmumps_c.h`` / ``xmumps_c.h`` (each carries
#      its own widened widths inline; a consumer ``#include <qmumps_c.h>``
#      needs only layers 1-2 for mumps_compat.h / mumps_int_def.h) plus
#      ``mumps_c_types_extended.h`` — the bridge-only force-include that
#      #includes the upstream mumps_c_types.h from layer 1 and overrides
#      its ``double`` widths. The extended header is shipped for
#      completeness/rebuilds; consumers of the wrappers do not include it.
if(NOT BASELINE_BUILD AND TARGET ${LIB_PREFIX}mumps)
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_include)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_include/
            DESTINATION include/mumps
            FILES_MATCHING PATTERN "*.h")
    endif()
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/c/include/mumps_int_def.h)
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/c/include/mumps_int_def.h
            DESTINATION include/mumps)
    endif()
    if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/target_${TARGET_NAME}/c/include)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/target_${TARGET_NAME}/c/include/
            DESTINATION include/mumps
            FILES_MATCHING PATTERN "*.h")
    endif()
endif()

if(NEEDS_MULTIFLOATS)
    # multifloats defines its own install rules; only install local helpers.
    # multifloats_mpi(_f) link MPI; tag them accordingly.
    set(_mf_mpi_helpers multifloats_mpi multifloats_mpi_f)
    foreach(_t la_constants_mw la_xisnan_mw multifloats_mpi multifloats_mpi_f)
        if(TARGET ${_t})
            set(_mf_mpi_arg "")
            if("${_t}" IN_LIST _mf_mpi_helpers)
                set(_mf_mpi_arg "MPI")
            endif()
            fortran_install_library(${_t} EXPORT ${_t}Targets NAMESPACE eplinalg:: ${_mf_mpi_arg})
            get_target_property(_fc_moddir ${_t} Fortran_MODULE_DIRECTORY)
            if(_fc_moddir)
                fortran_install_modules(${_t})
            endif()
        endif()
    endforeach()
endif()

if(NEEDS_QUAD_MPI)
    # quad_mpi / quad_mpi_f both link (or stand in for) MPI; tag them.
    foreach(_t quad_mpi quad_mpi_f)
        if(TARGET ${_t})
            fortran_install_library(${_t} EXPORT ${_t}Targets NAMESPACE eplinalg:: MPI)
            get_target_property(_fc_moddir ${_t} Fortran_MODULE_DIRECTORY)
            if(_fc_moddir)
                fortran_install_modules(${_t})
            endif()
        endif()
    endforeach()
endif()

# ── Tests ────────────────────────────────────────────────────────────
# Differential precision test suite. Each library has its own subdir
# under tests/ which is added when its tests/<lib>/CMakeLists.txt is
# present in the staging copy.
#
# We deliberately do NOT gate on BUILD_TESTING: the multifloats
# sub-build sets it locally to OFF to suppress its own tests, and
# that local value shadows the cache value from include(CTest) at
# the parent scope. Always-on with a directory check gives us the
# single signal "is the staging dir self-contained?" — drop tests/
# from staging to skip them.
option(MIGRATOR_BUILD_TESTS
    "Build the differential-precision test suite (tests/)" ON)
if(MIGRATOR_BUILD_TESTS AND IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests)
    enable_testing()
    add_subdirectory(tests)
endif()
