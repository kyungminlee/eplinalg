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
# with LIB_PAIR_PREFIX="" so the test framework links the standard archives
# (blas / lapack / pblas / scalapack) directly. We skip every
# add_migrated_* call and its post-link wiring under that mode — the
# bare archive names would otherwise collide with the standard archive
# of the same name (e.g. ``target_link_libraries(blas PUBLIC blas)``).
if("${LIB_PAIR_PREFIX}" STREQUAL "")
    set(BASELINE_BUILD TRUE)
else()
    set(BASELINE_BUILD FALSE)
endif()

# Release version: staging copies the repo-root VERSION file (single
# source of truth) into the staged tree; fall back for hand-rolled trees.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION")
    file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" EPLINALG_VERSION LIMIT_COUNT 1)
else()
    set(EPLINALG_VERSION 0.0.0)
endif()

project(${TARGET_NAME}_libraries VERSION ${EPLINALG_VERSION} LANGUAGES Fortran C CXX)

# ── Compiler-aware install infrastructure ────────────────────────────
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
include(FortranCompiler)

# ── Global compiler flags ────────────────────────────────────────────
# Compile every object -fPIC so the static archives can be repackaged
# as shared libraries (--whole-archive). Default-PIE toolchains emit
# -fPIE objects whose direct PC32 references to global data (Fortran
# module variables, BLACS/MPI-op state) are rejected when linking a
# .so. Negligible cost on x86-64; .a consumers are unaffected.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
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
endif()
if(MPI_TAG)
    message(STATUS "MPI tag: ${MPI_TAG}")
elseif(MPI_C_FOUND)
    message(STATUS "MPI found but could not identify flavor/version — MPI-dependent libraries will not be tag-suffixed")
endif()

# ── libmpiseq compatibility gate ─────────────────────────────────────
# libmpiseq's stubs (runtime/mpiseq/mpiseq_c_stubs.c et al.) are written against
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

# ── Sequential (libmpiseq) release tagging ───────────────────────────
# A "seq" deliverable ships the MUMPS stack linked against the in-tree
# libmpiseq stub rather than a real MPI. Its MPI-dependent archives must
# not collide with the real-MPI (``-${MPI_TAG}``) archives in a shared
# install prefix, so when MUMPS_LIBSEQ_RELEASE is ON those libraries are
# stamped with the flavor tag ``seq`` in place of the detected MPI_TAG.
# PT-Scotch is a distributed (real-MPI) ordering and is forced OFF in
# this column (see the MUMPS_HAVE_PTSCOTCH gate).
#
# libmpiseq's stubs are Intel-MPI-ABI-specific (see the gate above), so a
# seq release is only well-formed on an Intel-MPI toolchain; requesting it
# elsewhere is a hard error rather than a silent non-seq build.
option(MUMPS_LIBSEQ_RELEASE
    "Build the sequential libmpiseq release: tag MPI-dependent archives '-seq' \
instead of '-${MPI_TAG}', force PT-Scotch off, and install libmpiseq."
    OFF)
if(MUMPS_LIBSEQ_RELEASE AND NOT MPISEQ_IS_INTEL_MPI)
    message(FATAL_ERROR
        "MUMPS_LIBSEQ_RELEASE requires an Intel-MPI toolchain (libmpiseq's "
        "stub ABI). Configure with the linux-impi preset, or leave the option "
        "OFF for a real-MPI (-${MPI_TAG}) build.")
endif()

# MPI_LIB_TAG — the flavor tag stamped onto MPI-dependent library
# filenames and export files. Normally the detected MPI flavor; ``seq``
# for the libmpiseq release. Deliberately distinct from the raw MPI_TAG,
# which must keep reporting the true MPI flavor for the mpiseq gate above
# and external-ScaLAPACK resolution below.
if(MUMPS_LIBSEQ_RELEASE)
    set(MPI_LIB_TAG "seq")
else()
    set(MPI_LIB_TAG "${MPI_TAG}")
endif()
if(MPI_LIB_TAG AND NOT "${MPI_LIB_TAG}" STREQUAL "${MPI_TAG}")
    message(STATUS "MPI library tag: ${MPI_LIB_TAG} (overriding ${MPI_TAG})")
endif()

# ── Multifloats (conditional) ────────────────────────────────────────
# _multifloats_use_binary_release() — default acquisition path. Downloads the
# upstream multifloats binary release and wires up the same two targets the
# source build produced, so nothing downstream changes:
#   * multifloats::multifloats  — imported C++ core, pinned to the portable
#     non-LTO archive (the only core whose weak/COMDAT template ops fold
#     against our non-LTO -O2 m/w archives).
#   * multifloatsf              — Fortran module, compiled HERE from the
#     release's shipped, pre-expanded source (a .mod is version-locked, so it
#     must be produced by our exact Fortran toolchain), keeping the
#     eplinalg::multifloatsf install/export path identical to the source build.
# A macro (not a function) so the imported/normal targets and the
# MULTIFLOATS_BINARY_PREFIX signal land in this directory scope directly.
#
# Knobs (all optional):
#   MULTIFLOATS_VERSION  release tag, default the pinned latest; "latest"
#                        resolves live via the GitHub API (falls back to pin).
#   MULTIFLOATS_ARCHIVE  local path / custom URL to a prefetched bundle
#                        (offline builds); overrides the GitHub asset URL.
macro(_multifloats_use_binary_release)
    include(FetchContent)
    include(CheckFortranSourceCompiles)

    # 1. Resolve the release version.
    set(_mf_pinned "v0.8.2")   # latest known good — bump as releases land
    if(NOT DEFINED MULTIFLOATS_VERSION OR MULTIFLOATS_VERSION STREQUAL "")
        set(MULTIFLOATS_VERSION "${_mf_pinned}"
            CACHE STRING "multifloats binary-release tag ('latest' to auto-resolve)")
    endif()
    if(MULTIFLOATS_VERSION STREQUAL "latest")
        set(_mf_json "${CMAKE_CURRENT_BINARY_DIR}/_mf_latest.json")
        file(DOWNLOAD
            "https://api.github.com/repos/kyungminlee/multifloats/releases/latest"
            "${_mf_json}" STATUS _mf_dl_status TIMEOUT 20)
        list(GET _mf_dl_status 0 _mf_dl_rc)
        set(_mf_resolved "")
        if(_mf_dl_rc EQUAL 0)
            file(READ "${_mf_json}" _mf_json_body)
            string(JSON _mf_resolved ERROR_VARIABLE _mf_json_err
                GET "${_mf_json_body}" tag_name)
        endif()
        if(_mf_resolved AND NOT _mf_resolved STREQUAL "tag_name-NOTFOUND")
            set(MULTIFLOATS_VERSION "${_mf_resolved}")
            message(STATUS "multifloats: resolved latest release = ${MULTIFLOATS_VERSION}")
        else()
            set(MULTIFLOATS_VERSION "${_mf_pinned}")
            message(WARNING "multifloats: could not resolve 'latest' "
                "(offline?); falling back to pin ${_mf_pinned}")
        endif()
    endif()

    # 2. Resolve the platform asset (linux-x86_64 | darwin-arm64).
    string(TOLOWER "${CMAKE_SYSTEM_NAME}" _mf_os)
    set(_mf_arch "${CMAKE_SYSTEM_PROCESSOR}")
    if(_mf_arch MATCHES "^(x86_64|amd64|AMD64)$")
        set(_mf_arch "x86_64")
    elseif(_mf_arch MATCHES "^(arm64|aarch64|ARM64)$")
        set(_mf_arch "arm64")
    endif()
    set(_mf_asset "libmultifloats-${MULTIFLOATS_VERSION}-${_mf_os}-${_mf_arch}.tar.gz")

    # 3. Locate the archive (local override or GitHub asset URL).
    if(DEFINED MULTIFLOATS_ARCHIVE AND NOT MULTIFLOATS_ARCHIVE STREQUAL "")
        set(_mf_url "${MULTIFLOATS_ARCHIVE}")
    else()
        set(_mf_url "https://github.com/kyungminlee/multifloats/releases/download/${MULTIFLOATS_VERSION}/${_mf_asset}")
    endif()
    message(STATUS "multifloats: binary release ${MULTIFLOATS_VERSION} — ${_mf_url}")

    # 4. Download + extract into the build tree. DOWNLOAD_EXTRACT_TIMESTAMP
    #    stamps extracted files with the archive's own mtimes (CMP0135) — it
    #    silences the dev warning and is harmless here since we key off content,
    #    not timestamps.
    FetchContent_Declare(multifloats_bin URL "${_mf_url}"
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_Populate(multifloats_bin)
    set(_mf_prefix "${multifloats_bin_SOURCE_DIR}")

    # 5. C++ core: force the portable non-LTO build. The core Config prefers an
    #    LTO-tagged archive when the consumer compiler matches, but LTO-IR
    #    COMDAT sections can't fold against our non-LTO m/w objects → include
    #    the nolto targets file directly to pin libmultifloats-nolto.a.
    set(_mf_core_nolto
        "${_mf_prefix}/lib/cmake/multifloats/multifloatsTargets-nolto.cmake")
    if(NOT EXISTS "${_mf_core_nolto}")
        message(FATAL_ERROR
            "multifloats binary release ${MULTIFLOATS_VERSION} ships no portable "
            "non-LTO core (${_mf_core_nolto}); its LTO cores cannot link against "
            "our non-LTO m/w archives. Pick another MULTIFLOATS_VERSION or build "
            "from source (-DMULTIFLOATS_FROM_SOURCE=ON).")
    endif()
    include("${_mf_core_nolto}")

    # Existing eplinalg targets link the header-only C++ core by its BARE name
    # via $<BUILD_INTERFACE:multifloats> (e.g. runtime/multifloats-mpi's bridge
    # includes <multifloats.h>). The source build creates that bare
    # `multifloats` target directly; the release ships only the namespaced
    # imported `multifloats::multifloats`. Provide a bare INTERFACE target
    # forwarding to it so those references resolve (include dirs + link)
    # unchanged. BUILD_INTERFACE-only usage → never leaks into an install export.
    if(NOT TARGET multifloats)
        add_library(multifloats INTERFACE)
        target_link_libraries(multifloats INTERFACE multifloats::multifloats)
    endif()

    # 6. Fortran module: compile from the shipped, pre-expanded source. Mirror
    #    the release Config's own source-build branch (quad/noquad by REAL(16)
    #    support; strict IEEE / no FMA contraction — the DD error-free
    #    transforms require it).
    check_fortran_source_compiles("
        program t
          real(16) :: x
          x = 1.0_16
        end program
    " _MF_HAS_REAL16 SRC_EXT f90)
    if(_MF_HAS_REAL16)
        set(_mf_src "${_mf_prefix}/share/multifloatsf/src/multifloats-quad.f90")
    else()
        set(_mf_src "${_mf_prefix}/share/multifloatsf/src/multifloats-noquad.f90")
    endif()
    if(NOT EXISTS "${_mf_src}")
        message(FATAL_ERROR
            "multifloats binary release: Fortran source not found at ${_mf_src} "
            "(the bundle may be incomplete).")
    endif()

    add_library(multifloatsf STATIC "${_mf_src}")
    set_target_properties(multifloatsf PROPERTIES POSITION_INDEPENDENT_CODE ON)
    # Private module dir (NOT the shared fmod tree): fortran_install_modules
    # installs the target's whole Fortran_MODULE_DIRECTORY, so a private dir
    # keeps the multifloatsf package to just multifloats.mod. INSTALL_INTERFACE
    # points at the shared install location every eplinalg Fortran package uses.
    set(_mf_moddir "${PROJECT_BINARY_DIR}/_mf_fmod/${FORTRAN_MOD_COMPAT_TAG}")
    file(MAKE_DIRECTORY "${_mf_moddir}")
    set_target_properties(multifloatsf PROPERTIES
        Fortran_MODULE_DIRECTORY "${_mf_moddir}")
    _fc_make_relative_installdir(_mf_rel_libdir "${CMAKE_INSTALL_LIBDIR}")
    target_include_directories(multifloatsf PUBLIC
        $<BUILD_INTERFACE:${_mf_moddir}>
        $<INSTALL_INTERFACE:${_mf_rel_libdir}/fmod/${FORTRAN_MOD_COMPAT_TAG}>)
    target_link_libraries(multifloatsf PUBLIC multifloats::multifloats)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
        target_compile_options(multifloatsf PRIVATE -fp-model=precise)
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|LLVMFlang|Flang")
        target_compile_options(multifloatsf PRIVATE -ffp-contract=off)
    endif()

    # Signal EplinalgInstall to install a nolto-only passthrough of the core
    # (the sub-build used to self-install it; the binary release does not).
    set(MULTIFLOATS_BINARY_PREFIX "${_mf_prefix}" CACHE INTERNAL
        "Extracted multifloats binary-release prefix (core passthrough install)")
endmacro()

if(NEEDS_MULTIFLOATS)
    # ── How the multifloats double-double dependency is obtained ─────────
    # Three modes, in priority order:
    #   1. MULTIFLOATS_DIR defined      → local source checkout (add_subdirectory)
    #   2. MULTIFLOATS_FROM_SOURCE=ON   → git FetchContent source build (legacy)
    #   3. (default) binary release     → download the upstream prebuilt bundle
    #
    # The binary release is the DEFAULT. It ships the C++ double-double core as
    # the portable, untagged non-LTO libmultifloats-nolto.a — the only core
    # build whose weak/COMDAT template ops (float64x2 operators, deduped into
    # every m/w archive we compile non-LTO -O2) fold cleanly at link time. An
    # LTO core archive carries LTO-IR COMDAT sections that ld cannot fold
    # against our native objects → "defined in discarded section". We consume
    # that core as a binary and compile the Fortran module `multifloatsf`
    # ourselves from the release's shipped, pre-expanded source, so its
    # version-locked .mod stays under our exact toolchain and the
    # eplinalg::multifloatsf install/export path is unchanged.
    if(DEFINED MULTIFLOATS_DIR OR MULTIFLOATS_FROM_SOURCE)
        find_program(FYPP_EXECUTABLE fypp REQUIRED)
        # Use the multifloats top-level directory (not src/ and fsrc/
        # separately): upstream's src/CMakeLists.txt uses CMAKE_SOURCE_DIR
        # which is only correct when multifloats is the top-level project.
        #
        # We deliberately do NOT pass EXCLUDE_FROM_ALL. The C++ core
        # ``multifloats`` is a genuine runtime dependency of the m/w stack (its
        # extern "C" adddd / cadddd / … are called by the Fortran module
        # procedures and are not baked into any archive), so a find_package
        # consumer must be able to link it and resolve ``multifloats::multifloats``
        # (referenced by multifloatsfTargets). The sub-build already installs it
        # as package ``multifloats`` (namespace ``multifloats::``);
        # EXCLUDE_FROM_ALL suppresses that install *action* while still
        # registering the export set, which both hides the archive from the
        # prefix and, if we add our own, trips CMake's "target in multiple
        # export sets" error. Letting the sub-build install is the clean path.
        # Its own multifloatsf install stays gated off
        # (MULTIFLOATSF_INSTALL_PRECOMPILED_MOD), so we still install
        # multifloatsf ourselves below with no conflict.
        #
        # BUILD_TESTING/MULTIFLOATS_BUILD_BENCH are set as plain (non-cache)
        # variables so they only suppress tests/benches *inside* the multifloats
        # sub-build (leaving just the two libraries) — without EXCLUDE_FROM_ALL
        # these gates keep the sub-build lean. A `CACHE … FORCE` here would also
        # disable our own `tests/` directory at the end of this file.
        set(BUILD_TESTING OFF)
        set(MULTIFLOATS_BUILD_BENCH OFF)
        # Build the C++ core WITHOUT LTO (see the COMDAT note above; a non-LTO
        # build installs the portable, untagged libmultifloats.a so any consumer
        # resolves multifloats::multifloats regardless of its toolchain). Plain
        # (non-cache) set() → CMP0077 lets it win over the sub-build's option().
        set(MULTIFLOATS_USE_LTO OFF)
        if(DEFINED MULTIFLOATS_DIR)
            message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
            add_subdirectory(${MULTIFLOATS_DIR}
                ${CMAKE_CURRENT_BINARY_DIR}/_mf)
        else()
            include(FetchContent)
            set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
                CACHE STRING "Git URL for the multifloats library")
            # v0.6.0 is the floor: earlier tags leak ${CMAKE_SOURCE_DIR}
            # include paths when add_subdirectory'd (upstream issue #23) and
            # are not supported here.
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
                ${CMAKE_CURRENT_BINARY_DIR}/_mf)
        endif()
    else()
        _multifloats_use_binary_release()
    endif()

    # multifloats_mpi / multifloats_mpi_f — first-party MPI bridge for the
    # real64x2 arithmetic, a standalone library under runtime/. Its own
    # CMakeLists.txt defines both targets (STATIC when MPI_CXX_FOUND, else
    # an INTERFACE header-only target); it relies on this parent scope
    # having already found MPI and created the `multifloats` target.
    # staging.py copies runtime/multifloats-mpi/ into the staged tree, so
    # the relative path resolves there.
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/runtime/multifloats-mpi/CMakeLists.txt")
        add_subdirectory(runtime/multifloats-mpi)
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
    # quad_mpi / quad_mpi_f — first-party MPI bridge for the kind16
    # (__float128) arithmetic, a standalone library under runtime/ with its
    # own CMakeLists.txt defining both targets. staging.py copies
    # runtime/quad-mpi/ into the staged tree, so the relative path
    # resolves there.
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/runtime/quad-mpi/CMakeLists.txt")
        add_subdirectory(runtime/quad-mpi)
    endif()
endif()

# ── Helper functions ─────────────────────────────────────────────────
# add_standard_*/add_migrated_* library assembly plus the small shared
# helpers (fortran_relax_line_length, compile_c_as_cxx, link_if_present).
# The privatize gate must be defined before the helpers call it from
# add_migrated_* for PRIVATIZED_LIBRARIES entries.
include(${CMAKE_CURRENT_SOURCE_DIR}/EplinalgPrivatizeGate.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/EplinalgLibraryHelpers.cmake)

# ── Libraries in dependency order ────────────────────────────────────

# 1. BLAS (no deps).
#
# Standard-precision archive ``blas`` is built from the vendored Netlib
# BLAS sources (staged into ``_refblas_src/`` by ``migrator stage``).
# The migrated archive ``${LIB_PAIR_PREFIX}blas`` carries only Q/X/E/Y/M/W
# clones and PUBLIC-links the standard archive so downstreams resolve
# the S/D/C/Z entry points (DGEMM_, ZGEMM_, ...) through one link line.
# Cross-precision routines (dsdot, sdsdot) are excluded — the migrator
# also skips them (see codegen/recipes/blas.yaml).
add_standard_fortran_library(blas
    ${CMAKE_CURRENT_SOURCE_DIR}/_refblas_src
    EXCLUDE_REGEX "/(dsdot|sdsdot)\\.f9?0?$")
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(blas)
    link_if_present(${LIB_PAIR_PREFIX}blas blas)
endif()

# 1.5. XBLAS (depends on BLAS).
#
# Netlib's Extra-Precise BLAS (1.0.248) is a C library that exposes
# Level 1/2/3 BLAS-like entry points with an additional ``prec``
# parameter selecting the internal accumulation precision
# (single / double / extra). The migrated archive ``${LIB_PAIR_PREFIX}xblas``
# carries Q/X/E/Y clones (BLAS_qgemv_x, BLAS_xgemv_x, ...) of every
# uniform-precision _x routine; mixed-input variants
# (BLAS_dgemv_d_s_x, ...) are skipped by the recipe — they're not
# expressible in a single extended-precision target and unreachable
# from LAPACK.
#
# There is no standard-precision archive for XBLAS in this build:
# nothing in the project links the regular C XBLAS (LAPACK's
# *_extended.f routines that would call into it are not built).
if(NOT BASELINE_BUILD)
    add_migrated_c_library(xblas)
    if(TARGET ${LIB_PAIR_PREFIX}xblas)
        # f2c-bridge.h defaults to no Fortran-symbol mangling. gfortran
        # emits lowercase + trailing underscore (the "Add_" convention),
        # so define CONFIG_FC_UNDERSCORE to make the bridge produce
        # ``blas_qgemv_x_`` etc. — matching what Fortran callers expect.
        target_compile_definitions(${LIB_PAIR_PREFIX}xblas PRIVATE
            CONFIG_FC_UNDERSCORE)
        if(TARGET xblas_common)
            target_compile_definitions(xblas_common PRIVATE
                CONFIG_FC_UNDERSCORE)
        endif()
        link_if_present(${LIB_PAIR_PREFIX}xblas ${LIB_PAIR_PREFIX}blas)
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
    # multifloats_mpi / quad_mpi only exist under their respective
    # NEEDS_* modes, so the TARGET checks inside link_if_present gate
    # the bridge links to the matching arithmetic.
    link_if_present(${LIB_PAIR_PREFIX}blacs multifloats_mpi quad_mpi blacs)
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
    # The xx-family iterative-refinement helpers (qla_gerfsx_extended,
    # xla_gerfsx_extended, ...) live inside ${LIB_PAIR_PREFIX}lapack and call
    # XBLAS extra-precise Level-2 entry points (BLAS_QGEMV_X_, ...). With
    # static archives, linking ${LIB_PAIR_PREFIX}lapack into a test executable
    # pulls those *_extended objects via undefined-symbol resolution as
    # soon as any xx driver (qgesvxx / xgesvxx / qgbsvxx / ...) is reached
    # transitively, so the XBLAS symbols must be resolvable on the same
    # link line — the xblas dependency is PUBLIC so every consumer of
    # the migrated lapack archive picks it up automatically.
    link_if_present(${LIB_PAIR_PREFIX}lapack
        ${LIB_PAIR_PREFIX}blas lapack ${LIB_PAIR_PREFIX}xblas)
    # The per-target precision helpers (la_constants_{ey,qx,mw},
    # la_xisnan_{ey,qx,mw}) ride in the PREFIXED ${LIB_PAIR_PREFIX}lapack
    # archive (routed there by staging.py), not the shared lapack_common.
    # They ``USE`` the standard la_constants / la_xisnan modules supplied
    # by the std ``lapack`` archive; the ${LIB_PAIR_PREFIX}lapack -> lapack
    # PUBLIC link above already resolves ``la_xisnan.mod`` at compile
    # time. The lapack_common -> lapack link below is defensive — it
    # only fires if some other precision-independent member keeps a
    # lapack_common archive alive.
    link_if_present(lapack_common lapack)
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
    link_if_present(${LIB_PAIR_PREFIX}ptzblas ${LIB_PAIR_PREFIX}blas ptzblas)
endif()

# 5. PBBLAS (depends on BLAS, BLACS; scalapack_common wired post-section-8).
#
# Standard archive ``pbblas`` is built from PBLAS/PBBLAS/. Pure Fortran;
# no extra dependencies beyond BLAS+BLACS at link time.
add_standard_fortran_library(pbblas
    ${CMAKE_CURRENT_SOURCE_DIR}/_pbblas_src)
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(pbblas)
    # PBBLAS also calls NUMROC/ICEIL/ILCM from scalapack_common, which
    # is created later (section 7's add_migrated_fortran_library), so
    # that PUBLIC link is wired after section 8 — see the "PBBLAS →
    # scalapack_common" block below.
    link_if_present(${LIB_PAIR_PREFIX}pbblas
        ${LIB_PAIR_PREFIX}blas ${LIB_PAIR_PREFIX}blacs pbblas)
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
    link_if_present(${LIB_PAIR_PREFIX}pblas
        ${LIB_PAIR_PREFIX}blas ${LIB_PAIR_PREFIX}blacs
        ${LIB_PAIR_PREFIX}ptzblas ${LIB_PAIR_PREFIX}pbblas pblas)
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
    # scalapack_common (which carries NUMROC/ICEIL/ILCM) is linked
    # automatically by add_migrated_fortran_library(scalapack) above.
    link_if_present(${LIB_PAIR_PREFIX}scalapack
        ${LIB_PAIR_PREFIX}lapack ${LIB_PAIR_PREFIX}blacs ${LIB_PAIR_PREFIX}pblas scalapack)
endif()

# 8. ScaLAPACK C entry points — folded into the scalapack archive.
#
# The C entry points under SRC/*.c (clamov, dlamov, pdlaiect, pcrot,
# ...) plus REDIST/SRC/*.c (pdgemr2d_, pztrmr2d_, ...) and TOOLS/*.c
# (sl_gridreshape, reshape) are mutually dependent with ScaLAPACK's
# Fortran routines (which call back into qlamov_, pqgemr2d_, ...). A
# separate C archive would make that a target cycle — which breaks
# ``find_package`` of the exported targets — so the C objects fold
# directly into the scalapack archive, keeping the reference INTERNAL
# to one archive: ``ld`` resolves intra-archive back-references via
# the archive symbol index regardless of member order (no
# ``--start-group`` needed) and a single ``eplinalg::scalapack``
# carries both halves. Built as OBJECT libraries to keep the C (and,
# under multifloats, C-as-C++) compile settings off the Fortran
# sources. ``-DAdd_`` selects the lowercase + trailing-underscore
# Fortran symbol mangling that gfortran emits (upstream
# SLmake.inc.example's convention).
if(TARGET scalapack)
    file(GLOB _scalapack_c_std
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src/*.c)
    # SRC/pzaxpy.c, pzdotc.c, pzdotu.c are PBLAS-v1-era leftovers that
    # upstream ScaLAPACK never compiles (its SRC/CMakeLists.txt only builds
    # pcrot.c/pzrot.c from the C files). Modern ANSI C definitions of
    # pzaxpy_/pzdotc_/pzdotu_ live in PBLAS/SRC (libpblas); the SRC/ copies
    # are older K&R duplicates that reference pberror_/ptop, which exist
    # only in PBLAS TESTING/TIMING and would shadow the good definitions
    # from inside the same archive. Mirrors codegen/recipes/scalapack_c.yaml
    # skip_files.
    list(FILTER _scalapack_c_std EXCLUDE REGEX "/(pzaxpy|pzdotc|pzdotu)\\.c$")
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
    # scalapack carries C objects that call MPI; propagate the link
    # privately for its own objects and interface-only for downstream
    # Fortran.
    if(MPI_C_FOUND)
        target_link_libraries(scalapack PRIVATE MPI::MPI_C)
        target_link_libraries(scalapack INTERFACE $<LINK_ONLY:MPI::MPI_C>)
    endif()
endif()
if(NOT BASELINE_BUILD AND TARGET ${LIB_PAIR_PREFIX}scalapack)
    # Migrated C clones (qlamov_, pqgemr2d_, ...) fold into the migrated
    # scalapack archive pair instead of a separate qscalapack_c. Read the
    # scalapack_c manifest for the precision/common C split and build each
    # as an OBJECT lib; the precision objects fold into
    # ${LIB_PAIR_PREFIX}scalapack and the arithmetic-independent common object
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
                    compile_c_as_cxx(${_obj})
                    if(TARGET multifloats_mpi)
                        target_link_libraries(${_obj} PRIVATE multifloats_mpi)
                    endif()
                endif()
            endforeach()
        endif()
        if(TARGET scalapack_c_prec_obj)
            target_sources(${LIB_PAIR_PREFIX}scalapack
                PRIVATE $<TARGET_OBJECTS:scalapack_c_prec_obj>)
        endif()
        if(TARGET scalapack_c_common_obj)
            if(TARGET scalapack_common)
                target_sources(scalapack_common
                    PRIVATE $<TARGET_OBJECTS:scalapack_c_common_obj>)
            else()
                target_sources(${LIB_PAIR_PREFIX}scalapack
                    PRIVATE $<TARGET_OBJECTS:scalapack_c_common_obj>)
            endif()
        endif()
        # The migrated scalapack archive carries the C objects: expose
        # the C source dir (internal prototype headers) and propagate
        # MPI + the multifloats bridge.
        target_include_directories(${LIB_PAIR_PREFIX}scalapack
            PUBLIC $<BUILD_INTERFACE:${_scalapack_c_dir}/src>)
        if(MPI_C_FOUND)
            target_link_libraries(${LIB_PAIR_PREFIX}scalapack PRIVATE MPI::MPI_C)
            target_link_libraries(${LIB_PAIR_PREFIX}scalapack
                INTERFACE $<LINK_ONLY:MPI::MPI_C>)
        endif()
        if(C_AS_CXX)
            if(TARGET multifloats)
                # Header-only C++ core — build-time only (symbols baked into
                # the archive). Keep it out of the install interface so
                # find_package(mwscalapack) doesn't chase an uninstalled core.
                target_link_libraries(${LIB_PAIR_PREFIX}scalapack PUBLIC $<BUILD_INTERFACE:multifloats>)
            endif()
            if(TARGET multifloats_mpi)
                target_link_libraries(${LIB_PAIR_PREFIX}scalapack PUBLIC multifloats_mpi)
            endif()
        endif()
    endif()
endif()

# PBBLAS → scalapack_common (NUMROC / ICEIL / ILCM).
#
# pb*tran / pb*trnv / pb*trsrt / pb*trget / pb*trst1 in PBBLAS call
# NUMROC / ICEIL / ILCM. Those three precision-independent integer
# helpers live in the family-independent scalapack_common (routed there
# by scalapack.yaml), never in a prefixed archive: a per-family copy
# would duplicate three symbols MKL also provides and multiply-define
# against MKL's numroc_/iceil_/ilcm_ under extended-MKL, whereas a
# consumer simply drops scalapack_common when MKL supplies that layer
# (same profile as blacs_common). scalapack_common is created in
# section 7 (after PBBLAS, section 5), so the PUBLIC link is wired here
# where the target is guaranteed to exist. tests/pbblas links only
# ${LIB_PAIR_PREFIX}pbblas — it transitively pulls scalapack_common for
# these symbols without dragging in the full ${LIB_PAIR_PREFIX}scalapack.
#
# Wrapped in $<BUILD_INTERFACE:> so the dependency serves the in-tree
# tests but is NOT written into install(EXPORT). scalapack_common lives
# in a DIFFERENT package (eyscalapack's export, per _install_library_pair)
# than eypbblas (pbblas' export); leaking it into eypbblas'
# INTERFACE_LINK_LIBRARIES makes find_package(eypbblas) reference an
# imported target defined only by find_package(eyscalapack) — and a
# consumer that finds eypbblas first then trips eypbblas_FOUND = FALSE.
# Keeping it build-only matches the repo convention (no cross-package
# find_dependency; consumers list every archive on their own link line,
# resolving NUMROC/ICEIL/ILCM from MKL or from eplinalg::scalapack_common
# via find_package(eyscalapack)) — identical to how the MUMPS parallel
# stack is BUILD_INTERFACE-wrapped out of the solver exports.
if(NOT BASELINE_BUILD AND TARGET ${LIB_PAIR_PREFIX}pbblas AND TARGET scalapack_common)
    target_link_libraries(${LIB_PAIR_PREFIX}pbblas PUBLIC $<BUILD_INTERFACE:scalapack_common>)
endif()

# 9. MUMPS (depends on LAPACK, BLAS, ScaLAPACK; needs MPI Fortran) —
# the largest section: orderings (PORD/METIS/Scotch/PT-Scotch), the
# typed q/x solvers, the C bridge, and the genuine dz/sc solvers.
include(${CMAKE_CURRENT_SOURCE_DIR}/EplinalgMumps.cmake)

# 10. libmpiseq — sequential MPI stub from MUMPS 5.8.2 upstream, a
# standalone library under runtime/ like the multifloats-mpi and
# quad-mpi bridges. Its CMakeLists (runtime/mpiseq/) holds the full
# ABI-matching story; the first-party stub sources live next to it.
# The upstream libseq sources are staged into ``_mpiseq_src/`` by
# ``migrator stage`` (see ``_std_dirs`` in ``codegen/migrator/staging.py``);
# the build silently skips when that directory is absent (e.g., a bare
# repo with no upstream MUMPS vendored). Intel-MPI-only — see the
# MPISEQ_IS_INTEL_MPI gate above.
if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/_mpiseq_src AND MPISEQ_IS_INTEL_MPI
        AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/runtime/mpiseq/CMakeLists.txt")
    add_subdirectory(runtime/mpiseq)
endif()

# ── Install ──────────────────────────────────────────────────────────
# Export sets, Config packages, and install rules for every archive
# built above.
include(${CMAKE_CURRENT_SOURCE_DIR}/EplinalgInstall.cmake)

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
