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

    # la_constants_mf helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_constants_mf.f90")
        add_library(la_constants_mf STATIC "${MF_HELPERS_DIR}/la_constants_mf.f90")
        fortran_module_layout(la_constants_mf)
        target_link_libraries(la_constants_mf PUBLIC $<BUILD_INTERFACE:multifloatsf>)
    endif()

    # la_xisnan_mf helper module (needs Fortran multifloats module)
    if(EXISTS "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        add_library(la_xisnan_mf STATIC "${MF_HELPERS_DIR}/la_xisnan_mf.f90")
        fortran_module_layout(la_xisnan_mf)
        target_link_libraries(la_xisnan_mf PUBLIC $<BUILD_INTERFACE:multifloatsf>)
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
# ``add_standard_fortran_library`` for C libraries (BLACS, PBLAS,
# ScaLAPACK_C). Includes MPI privately and propagates it via INTERFACE
# for downstream Fortran link, mirroring add_migrated_c_library.
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
        if(TARGET la_constants_mf)
            target_link_libraries(${_precision_target} PUBLIC la_constants_mf)
        endif()
        if(TARGET la_xisnan_mf)
            target_link_libraries(${_precision_target} PUBLIC la_xisnan_mf)
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

    # Dual-interface compile: files listed in
    # ``${lib_name}_DUAL_INTERFACE_SOURCES`` are already compiled once
    # above (they live in _precision / _common, since every such file
    # is also one of BLACS's regular sources) and expose the
    # Fortran-callable entry points by default. Compile them a second
    # time as an OBJECT library with ``-DCallFromC`` so the C-callable
    # entry points (``Cblacs_*``, ``Cdgesd2d`` etc.) are also present
    # in the final static archive. Empty ``_dual`` list is a no-op.
    if(_dual)
        set(_dual_target "${_precision_target}_c_api")
        add_library(${_dual_target} OBJECT ${_dual})
        target_compile_definitions(${_dual_target} PRIVATE CallFromC)
        target_include_directories(${_dual_target} PRIVATE $<BUILD_INTERFACE:${_dir}/src>)
        if(MPI_C_FOUND)
            target_link_libraries(${_dual_target} PRIVATE MPI::MPI_C)
        endif()
        if(C_AS_CXX)
            set_source_files_properties(${_dual} PROPERTIES LANGUAGE CXX)
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
        target_sources(${_precision_target}
                       PRIVATE $<TARGET_OBJECTS:${_dual_target}>)
    endif()

    # Multifloats: compile C as C++ for operator overloading support.
    # ``copy_files``-staged Fortran helpers (e.g. pblas/pilaenv.f) ride in
    # COMMON_SOURCES but must retain their native language; filter them
    # out of the language override.
    #
    # ScaLAPACK_C is exempt: REDIST/SRC sources (pctrmr2.c, pddgemr2.c,
    # …) use multi-line K&R-style function definitions with embedded
    # block comments inside the parameter declaration list — patterns
    # the migrator's K&R-to-ANSI converter doesn't unwind. Plain C is
    # legal here: REDIST routines only do block-cyclic copies and
    # ``sizeof``-based pointer arithmetic, never per-element math on
    # multifloats values, so they don't need operator overloading.
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
# scalapack_c clones are produced by the migrated archive's
# DUAL_INTERFACE_SOURCES path — the std archive doesn't replicate them.
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
# helper modules (la_constants_ep / la_xisnan_ep / la_constants_mf /
# la_xisnan_mf) are excluded — they only make sense in the migrated
# archive since they reference quad / multifloats kinds.
add_standard_fortran_library(lapack
    ${CMAKE_CURRENT_SOURCE_DIR}/_reflapack_src
    EXCLUDE_REGEX "/(dsgesv|zcgesv|dsposv|zcposv|dsgels|zcgels|dlag2s|slag2d|zlag2c|clag2z|dlat2s|zlat2c|la_constants_ep|la_constants_mf|la_xisnan_ep|la_xisnan_mf)\\.([fF]|[fF]90)$")
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
    # The migrated lapack_common archive holds the extended-precision-only
    # helpers (la_constants_ep, la_xisnan_ep) which ``USE`` the standard
    # la_constants / la_xisnan modules supplied by the std ``lapack``
    # archive. The link is needed so gfortran resolves ``la_xisnan.mod``
    # during compilation, not just at final link time.
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

# 4b. ScaLAPACK tools (NUMROC / ICEIL / ILCM).
#
# Three precision-independent integer helpers carved out of ScaLAPACK's
# TOOLS/ tree into their own tiny library. PBBLAS' pb*tran / pb*trnv /
# pb*trsrt / pb*trget / pb*trst1 reference NUMROC / ICEIL / ILCM, but
# tests/pbblas links only ${LIB_PREFIX}pbblas / blacs / blas — no
# scalapack — so the symbols had to come from somewhere reachable
# without dragging in the full ScaLAPACK archive.
#
# No standard counterpart archive: the baseline ``scalapack`` archive
# already absorbs TOOLS/*.f directly (see section 7), and tests/pbblas
# self-skips on baseline (HAVE_REAL16 is required). The migrated
# library carries the symbols for both ``${LIB_PREFIX}pbblas`` and
# ``${LIB_PREFIX}scalapack`` (which both PUBLIC-link it below), with
# the canonical un-prefixed Fortran ABI names (numroc_, iceil_, ilcm_)
# since the symbols are precision-independent and pass the migrator
# pipeline unrenamed.
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(scalapack_tools)
endif()

# 5. PBBLAS (depends on BLAS, BLACS, scalapack_tools).
#
# Standard archive ``pbblas`` is built from PBLAS/PBBLAS/. Pure Fortran;
# no extra dependencies beyond BLAS+BLACS at link time.
add_standard_fortran_library(pbblas
    ${CMAKE_CURRENT_SOURCE_DIR}/_pbblas_src)
if(NOT BASELINE_BUILD)
    add_migrated_fortran_library(pbblas)
    if(TARGET ${LIB_PREFIX}pbblas)
        foreach(_dep ${LIB_PREFIX}blas ${LIB_PREFIX}blacs ${LIB_PREFIX}scalapack_tools)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}pbblas PUBLIC ${_dep})
            endif()
        endforeach()
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
        foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas
                     ${LIB_PREFIX}scalapack_tools)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}scalapack PUBLIC ${_dep})
            endif()
        endforeach()
    endif()
    if(TARGET ${LIB_PREFIX}scalapack AND TARGET scalapack)
        target_link_libraries(${LIB_PREFIX}scalapack PUBLIC scalapack)
    endif()
endif()

# 8. ScaLAPACK C (depends on LAPACK, BLACS, PBLAS, ScaLAPACK).
#
# Standard archive ``scalapack_c`` is built from the small number of
# C entry points under SRC/*.c (clamov, dlamov, pdlaiect, pcrot, ...)
# plus REDIST/SRC/*.c (pdgemr2d_, pztrmr2d_, ...) and TOOLS/*.c
# (sl_gridreshape, reshape). ``-DAdd_`` selects the Fortran symbol
# mangling matching gfortran's emit convention.
add_standard_c_library(scalapack_c
    ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src
    COMPILE_DEFINITIONS Add_)
if(TARGET scalapack_c)
    file(GLOB _scalapack_c_redist
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_redist_src/*.c)
    file(GLOB _scalapack_c_tools
         CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src/*.c)
    if(_scalapack_c_redist OR _scalapack_c_tools)
        target_sources(scalapack_c PRIVATE
            ${_scalapack_c_redist} ${_scalapack_c_tools})
    endif()
    target_include_directories(scalapack_c PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_redist_src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/_scalapack_tools_src>)
endif()
if(NOT BASELINE_BUILD)
    add_migrated_c_library(scalapack_c)
    if(TARGET ${LIB_PREFIX}scalapack_c)
        # REDIST/SRC sources select their Fortran-call entry-point naming
        # convention (trailing underscore vs UpCase vs no-change) at compile
        # time via ``Add_`` / ``UpCase`` / ``NoChange``. Upstream's
        # SLmake.inc.example sets ``-DAdd_`` to get the lowercase + trailing
        # underscore convention that gfortran emits, so we match.
        target_compile_definitions(${LIB_PREFIX}scalapack_c PRIVATE Add_)
        if(TARGET scalapack_c_common)
            target_compile_definitions(scalapack_c_common PRIVATE Add_)
        endif()
        foreach(_dep ${LIB_PREFIX}lapack ${LIB_PREFIX}blacs ${LIB_PREFIX}pblas ${LIB_PREFIX}scalapack)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PREFIX}scalapack_c PUBLIC ${_dep})
            endif()
        endforeach()
        if(TARGET scalapack_c)
            target_link_libraries(${LIB_PREFIX}scalapack_c PUBLIC scalapack_c)
        endif()
        # Reverse dep: ScaLAPACK's Fortran routines call qlamov_/xlamov_
        # (and the redistribution entry points pqgemr2d_/pxgemr2d_) which
        # live in ${LIB_PREFIX}scalapack_c. We used to wire this as a
        # PUBLIC link on ${LIB_PREFIX}scalapack so downstream callers
        # could link a single library — but that put qscalapack_c into
        # qscalapack's INTERFACE_LINK_LIBRARIES, creating a cycle with
        # qscalapack_c → qscalapack. The cycle blocks ``find_package``
        # because qscalapack's Targets file then fails the
        # ``_cmake_expected_target`` check before qscalapack_c is
        # loaded. Consumers calling Fortran ScaLAPACK routines that
        # internally jump to C-side helpers must now explicitly link
        # both ``${LIB_PREFIX}scalapack`` and ``${LIB_PREFIX}scalapack_c``.
    endif()
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
    endif()
    if(TARGET mumps_common AND MPI_Fortran_FOUND)
        target_link_libraries(mumps_common PRIVATE MPI::MPI_Fortran)
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
    # precision-aware ``mumps_c_types.h`` shadow) come from
    # ``tests/mumps/target_${TARGET_NAME}/c/include/`` — produced by
    # the migrator at stage time. They install AFTER the upstream
    # headers so the precision-aware shadow ``mumps_c_types.h``
    # shadows upstream's plain-``double`` version.
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
            -include mumps_c_types.h)

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
            -include mumps_c_types.h)

        # mumps_addr.c stays at upstream double precision (the migrated
        # Fortran's xtools.F / qtools.F UPDATE_PROGRESS pass DP per
        # keep-kind manifest — compiling with the shadow types would
        # corrupt the 8-byte caller buffer).
        add_library(${LIB_PREFIX}mumps_c_addr_obj OBJECT ${_mumps_c_src}/mumps_addr.c)
        target_include_directories(${LIB_PREFIX}mumps_c_addr_obj PRIVATE
            ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
        target_compile_definitions(${LIB_PREFIX}mumps_c_addr_obj PRIVATE Add_)

        # Combined STATIC archive. Bundles the per-arithmetic objects
        # with the type-agnostic upstream runtime into ONE library so
        # consumers link a single ``${LIB_PREFIX}mumps_c``.
        add_library(${LIB_PREFIX}mumps_c STATIC
            $<TARGET_OBJECTS:${LIB_PREFIX}mumps_c_obj>
            $<TARGET_OBJECTS:${LIB_PREFIX_COMPLEX}mumps_c_obj>
            $<TARGET_OBJECTS:${LIB_PREFIX}mumps_c_addr_obj>
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
        target_include_directories(${LIB_PREFIX}mumps_c
            PUBLIC
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_target}>
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_shared}>
                $<BUILD_INTERFACE:${_mumps_c_inc}>
                $<BUILD_INTERFACE:${_mumps_c_src}>
                $<INSTALL_INTERFACE:include/mumps>)
        target_compile_definitions(${LIB_PREFIX}mumps_c PRIVATE Add_)
        target_compile_options(${LIB_PREFIX}mumps_c PRIVATE
            -include mumps_c_types.h)
        if(MPI_C_FOUND)
            target_link_libraries(${LIB_PREFIX}mumps_c PRIVATE MPI::MPI_C)
            target_link_libraries(${LIB_PREFIX}mumps_c INTERFACE
                $<LINK_ONLY:MPI::MPI_C>)
        endif()
        # Link against the migrated Fortran archive at consumer time —
        # qmumps_c's bridge symbols call qmumps_f77_ etc.
        target_link_libraries(${LIB_PREFIX}mumps_c PUBLIC ${LIB_PREFIX}mumps)

        # Umbrella INTERFACE target. The qmumps Fortran archive and
        # the qmumps_c bridge are mutually recursive at link time:
        # qmumps.a references mumps_size_c_ (in qmumps_c.a's
        # mumps_addr.c) and qmumps_c.a references qmumps_f77_ (in
        # qmumps.a). A plain left-to-right scan misses one or the
        # other. The ``${LIB_PREFIX}mumps_full`` target wraps both
        # archives in a ``LINK_GROUP:RESCAN`` so ld --start-group /
        # --end-group resolves the cycle automatically. Consumers
        # who use ``kind16_libraries::${LIB_PREFIX}mumps_full``
        # don't have to know about the static-archive ordering
        # subtlety. Requires CMake 3.24+ on the consumer side.
        add_library(${LIB_PREFIX}mumps_full INTERFACE)
        # Propagate include directories (mumps_full is a meta-target;
        # by itself it has no headers, so pull them from qmumps_c
        # which carries the per-target mumps_c.h shadow path).
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
            "$<LINK_GROUP:RESCAN,${PROJECT_NAME}::${LIB_PREFIX}mumps,${PROJECT_NAME}::${LIB_PREFIX}mumps_c>")
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
set(_MPI_DEPENDENT_LIBS blacs pblas pbblas scalapack scalapack_c mumps)


# install_library_headers(<lib_name>): install any ``*.h`` headers
# under ``${lib_name}/src/`` to ``include/${lib_name}/`` and extend
# every related target's INTERFACE include directories so consumers
# linking against the installed library pick up the include path via
# ``find_package``. Per-library subdirs are used because some header
# names collide across libraries (e.g. PBLAS and ScaLAPACK_C both
# ship a ``pblas.h`` with different contents).
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
# `blas` → `eplinalgStdBlas`, `scalapack_c` → `eplinalgStdScalapack_c`, etc.
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
    # la_constants_mf / la_xisnan_mf helper archives (see ~line 403),
    # each installed as its own Config package (EXPORT ${_t}Targets,
    # ~line 1439). Like eplinalgStdBlas they are transparent deps, so
    # consumers only need find_package(mblas) — without these the
    # generated Config references undefined imported targets and sets
    # mblas_FOUND = FALSE.
    if(NEEDS_MULTIFLOATS)
        foreach(_mf_helper la_constants_mf la_xisnan_mf)
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

    # MUMPS-special: bundle the C-side bridge (${LIB_PREFIX}mumps_c)
    # and the umbrella INTERFACE (${LIB_PREFIX}mumps_full) into the
    # same export so ``find_package(qmumps)`` brings the Fortran
    # archive, the C-callable entry-points (qmumps_c, xmumps_c),
    # and the umbrella that wraps them in a LINK_GROUP:RESCAN.
    if("${lib_name}" STREQUAL "mumps" AND TARGET ${LIB_PREFIX}mumps_c)
        fortran_install_library(${LIB_PREFIX}mumps_c EXPORT ${_export} NAMESPACE eplinalg:: MPI)
        install(TARGETS ${LIB_PREFIX}mumps_full
            EXPORT ${_export})
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

# MUMPS additionally ships C-API entry-point headers in three layers,
# installed in this order so the precision-aware shadow overrides
# upstream where both define the same filename:
#
#   1. ``_mumps_upstream_include/*.h`` — upstream's dmumps_c.h,
#      zmumps_c.h, dmumps_struc.h, mumps_compat.h, mumps_c_types.h
#      (the last is the ``double``-typed version, shadowed below).
#   2. ``tests/mumps/c/include/mumps_int_def.h`` — fixed 32-bit
#      MUMPS_INT stub (upstream's is generated from a .in template).
#   3. ``tests/mumps/target_${TARGET_NAME}/c/include/*.h`` — the
#      per-target wrappers ``qmumps_c.h`` / ``xmumps_c.h`` plus a
#      precision-aware ``mumps_c_types.h`` shadow with QREAL fields.
#      Installed LAST so the shadow overrides upstream's plain-double
#      ``mumps_c_types.h``. Consumers ``#include <qmumps_c.h>`` get
#      the right precision via the include order set up in that
#      wrapper.
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
    foreach(_t la_constants_mf la_xisnan_mf multifloats_mpi multifloats_mpi_f)
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
