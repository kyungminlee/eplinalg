"""CMake project generation for the single-library ``migrator build`` path.

Extracted verbatim from ``__main__.py`` (Cluster 6) as part of the migrator
file-restructuring refactor. Behaviour is unchanged; ``cmd_build`` imports
``_generate_cmake`` from here. Pure string templating — the only inputs are the
``TargetMode`` and the source file lists.
"""
import shutil
from pathlib import Path


def _generate_cmake(output_dir: Path, lib_name: str, target_mode,
                    common_files: list[str], precision_files: list[str],
                    language: str = 'fortran',
                    project_root: Path | None = None,
                    ref_sources: list[Path] | None = None):
    """Generate a self-contained CMakeLists.txt in the output directory."""
    pmap = target_mode.prefix_map
    real_pfx = pmap['R'].lower()
    precision_lib = f'{real_pfx}{lib_name}'
    common_lib = f'{lib_name}_common'

    common_list = '\n    '.join(sorted(common_files))
    precision_list = '\n    '.join(sorted(precision_files))
    ref_list = '\n    '.join(f'"{p}"' for p in (ref_sources or []))

    # Default path to the vendored Intel MPI headers. ``project_root``
    # is resolved at generation time, so the generated CMakeLists.txt
    # works when built from a fresh out-of-tree output directory.
    _impi_default = str(((project_root or Path.cwd())
                         / 'external' / 'impi-headers').resolve())

    if language == 'c':
        # When targeting multifloats, the migrated C sources `#include
        # "multifloats_bridge.h"`. Copy + patch the bridge header into
        # the output dir so the migrated sources find it on the include
        # path, mirroring what cmd_stage does for the shared driver.
        c_mf_link = ''
        c_mf_deps = ''
        if target_mode.module_name is not None:
            _root = project_root or Path.cwd()
            mf_local = _root / 'external' / 'multifloats-mpi'
            bridge_h_src = mf_local / 'multifloats_bridge.h'
            if bridge_h_src.is_file():
                helpers_dst = output_dir / '_helpers'
                helpers_dst.mkdir(exist_ok=True)
                staged = helpers_dst / bridge_h_src.name
                shutil.copy2(bridge_h_src, staged)
                # Guard the bridge header's `#include <mpi.h>` against
                # the C++ MPI bindings — without it, scalapack_c's
                # C-as-C++ build pulls thousands of mpicxx.h templates
                # into the migrator-injected ``extern "C" { … }`` wrap
                # and fails to link. Same patch cmd_stage applies.
                text = staged.read_text()
                if 'MPICH_SKIP_MPICXX' not in text:
                    text = text.replace(
                        '#include <mpi.h>',
                        '#define MPICH_SKIP_MPICXX 1\n'
                        '#define OMPI_SKIP_MPICXX 1\n'
                        '#include <mpi.h>',
                        1,
                    )
                    staged.write_text(text)
            c_mf_link = """
# multifloats: FetchContent (or local via -DMULTIFLOATS_DIR) so the
# migrated sources can link against ``libmultifloats.a`` (C++) and
# include ``multifloats_bridge.h`` (staged into ./_helpers/).
set(BUILD_TESTING OFF CACHE BOOL "Disable tests in fetched multifloats" FORCE)
set(MULTIFLOATS_BUILD_BENCH OFF CACHE BOOL "Disable benches in fetched multifloats" FORCE)
if(DEFINED MULTIFLOATS_DIR)
    message(STATUS "Using local multifloats: ${MULTIFLOATS_DIR}")
    add_subdirectory(${MULTIFLOATS_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
else()
    include(FetchContent)
    set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
        CACHE STRING "Git URL for the multifloats library")
    set(MULTIFLOATS_GIT_TAG "v0.6.0"
        CACHE STRING "Git tag/branch/commit for multifloats (>= v0.6.0)")
    FetchContent_Declare(multifloats_fetch
        GIT_REPOSITORY ${MULTIFLOATS_GIT_REPO}
        GIT_TAG        ${MULTIFLOATS_GIT_TAG}
    )
    FetchContent_Populate(multifloats_fetch)
    add_subdirectory(
        ${multifloats_fetch_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/_mf EXCLUDE_FROM_ALL)
endif()
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/_helpers)
"""
            c_mf_deps = f"""
if(TARGET multifloats)
    target_link_libraries({precision_lib} PUBLIC multifloats)
    if(TARGET {common_lib})
        target_link_libraries({common_lib} PUBLIC multifloats)
    endif()
endif()
# multifloats's bridge header (multifloats_bridge.h) uses C++ ``using``
# declarations to expose ``float64x2`` at file scope. Migrated .c
# bodies need a C++ translation unit for those typedefs to resolve.
# Reclassify the migrated sources to LANGUAGE CXX so they go through
# the C++ compiler while keeping their .c extension on disk.
set_source_files_properties(${{PRECISION_SOURCES}} PROPERTIES LANGUAGE CXX)
if(COMMON_SOURCES)
    set_source_files_properties(${{COMMON_SOURCES}} PROPERTIES LANGUAGE CXX)
endif()
# C-as-C++ flags (mirror add_migrated_c_library in the shared driver):
#   - cxx_std_17: multifloats requires it.
#   - -fpermissive: BLACS/PBLAS bodies have implicit ``void *`` →
#     ``char *`` casts (e.g. BI_iMPI_amn.c:12) that C accepts and C++
#     rejects. Tolerated on g++; clang++ does NOT honor this flag,
#     so a clang-based multifloats build still needs the explicit-cast
#     overrides in ``recipes/<lib>/mfc_overrides/``.
#   - -Wno-write-strings: silences `const char[]` -> `char *` literals.
#   - MPICH_SKIP_MPICXX / OMPI_SKIP_MPICXX: keep mpicxx.h's templates
#     out of the migrator-injected ``extern "C" {{ ... }}`` wrap.
target_compile_features({precision_lib} PRIVATE cxx_std_17)
target_compile_options({precision_lib} PRIVATE -fpermissive -Wno-write-strings)
target_compile_definitions({precision_lib} PRIVATE
    MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
if(TARGET {common_lib})
    target_compile_features({common_lib} PRIVATE cxx_std_17)
    target_compile_options({common_lib} PRIVATE -fpermissive -Wno-write-strings)
    target_compile_definitions({common_lib} PRIVATE
        MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
endif()
"""

        # Add CXX to project() languages when multifloats is in play —
        # multifloats's targets request cxx_std_17 features from the
        # embedding project. Harmless on KIND targets (small detect
        # cost, no sources compiled as C++).
        project_langs = 'C CXX' if target_mode.module_name is not None else 'C'

        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} {project_langs})

# --- Compiler flags ---
set(CMAKE_C_FLAGS "${{CMAKE_C_FLAGS}} -w")

# --- MPI: default to vendored Intel MPI headers ---
# ``external/impi-headers/`` ships mpi.h and mpif.h at the Intel MPI
# ABI so the build compiles against a stable MPI surface without
# requiring every contributor to install an MPI runtime. Link-time
# libraries still come from whichever MPI runtime the user provides
# (impi-rt / OpenMPI / MPICH — headers are ABI-compatible).
# Users who want a different MPI's *headers* can override IMPI_HEADERS.
if(NOT DEFINED IMPI_HEADERS)
    set(IMPI_HEADERS "{_impi_default}"
        CACHE PATH "Path to vendored Intel MPI headers")
endif()
include_directories(${{IMPI_HEADERS}})
find_package(MPI COMPONENTS C QUIET)
{c_mf_link}
# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

# Header include path
include_directories(${{CMAKE_CURRENT_SOURCE_DIR}}/src)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()
{c_mf_deps}
# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
"""
    else:
        # If multifloats, we need to link against the multifloats library
        # AND build the la_constants_mw / la_xisnan_mw helper modules
        # that the migrated source depends on for la_constants USE clauses.
        mf_link = ""
        mf_deps = ""
        if target_mode.module_name is not None:
            # Resolve absolute paths to external dependencies so the
            # generated CMakeLists.txt works from any output directory.
            # MF helpers (la_constants_mw.f90 / la_xisnan_mw.f90) live
            # under recipes/<lib>/mf_helpers/. Prefer the per-recipe
            # directory if present; otherwise fall back to the upstream
            # SRC dir which (historically) shipped the EP helpers.
            _root = project_root or Path.cwd()
            _per_recipe_mf = _root / 'recipes' / lib_name / 'mf_helpers'
            if _per_recipe_mf.is_dir():
                _helpers_default = str(_per_recipe_mf.resolve())
            else:
                _helpers_default = str((_root / 'external' / 'lapack-3.12.1' / 'SRC').resolve())
            # multifloats-mpi extras: Fortran-side MPI handle module
            # used by MUMPS (``USE multifloats_mpi_f``).
            _mf_mpi_dir = (_root / 'external' / 'multifloats-mpi').resolve()
            mf_link = f"""
# Fetch the multifloats library from GitHub (default) or use a local
# checkout via -DMULTIFLOATS_DIR=/path/to/multifloats. We add the
# multifloats *top-level* directory so its CMakeLists.txt runs — the
# src/CMakeLists.txt references ``CMAKE_SOURCE_DIR/include`` which is
# wrong outside a top-level build. Tests/benches are suppressed via
# cache variables set before the subdirectory add.
set(BUILD_TESTING OFF CACHE BOOL "Disable tests in fetched multifloats" FORCE)
set(MULTIFLOATS_BUILD_BENCH OFF CACHE BOOL "Disable benches in fetched multifloats" FORCE)
if(DEFINED MULTIFLOATS_DIR)
    message(STATUS "Using local multifloats: ${{MULTIFLOATS_DIR}}")
    add_subdirectory(${{MULTIFLOATS_DIR}}
        ${{CMAKE_CURRENT_BINARY_DIR}}/_mf EXCLUDE_FROM_ALL)
else()
    include(FetchContent)
    set(MULTIFLOATS_GIT_REPO "https://github.com/kyungminlee/multifloats.git"
        CACHE STRING "Git URL for the multifloats library")
    # multifloats v0.6.0 fixed the ${{CMAKE_SOURCE_DIR}} include-path
    # leak (upstream issue #23). Earlier tags fail at configure when
    # add_subdirectory'd; don't drop below this floor.
    set(MULTIFLOATS_GIT_TAG "v0.6.0"
        CACHE STRING "Git tag/branch/commit for multifloats (>= v0.6.0)")
    message(STATUS "Fetching multifloats from ${{MULTIFLOATS_GIT_REPO}} (${{MULTIFLOATS_GIT_TAG}})")
    FetchContent_Declare(multifloats_fetch
        GIT_REPOSITORY ${{MULTIFLOATS_GIT_REPO}}
        GIT_TAG        ${{MULTIFLOATS_GIT_TAG}}
    )
    FetchContent_Populate(multifloats_fetch)
    add_subdirectory(
        ${{multifloats_fetch_SOURCE_DIR}}
        ${{CMAKE_CURRENT_BINARY_DIR}}/_mf EXCLUDE_FROM_ALL)
endif()

# Build the la_constants_mw and la_xisnan_mw helper modules. These
# re-export multifloats's MF_* constants under the M/W-prefixed names
# that the migrated LAPACK source uses via its rewritten
# ``USE LA_CONSTANTS_MW`` clause.
set(MF_HELPERS_DIR "{_helpers_default}"
    CACHE PATH "Directory containing la_constants_mw.f90 / la_xisnan_mw.f90")
if(EXISTS "${{MF_HELPERS_DIR}}/la_constants_mw.f90")
    add_library(la_constants_mw STATIC
        "${{MF_HELPERS_DIR}}/la_constants_mw.f90")
    set_target_properties(la_constants_mw PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_constants_mw PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    # la_constants_mw.f90 does ``use multifloats, only: real64x2`` —
    # the Fortran module lives in ``multifloatsf`` (the Fortran half),
    # not the C++ ``multifloats`` target.
    if(TARGET multifloatsf)
        target_link_libraries(la_constants_mw PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
    endif()
    if(TARGET multifloats)
        target_link_libraries(la_constants_mw PUBLIC multifloats)
    endif()
endif()
if(EXISTS "${{MF_HELPERS_DIR}}/la_xisnan_mw.f90")
    add_library(la_xisnan_mw STATIC
        "${{MF_HELPERS_DIR}}/la_xisnan_mw.f90")
    set_target_properties(la_xisnan_mw PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(la_xisnan_mw PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(TARGET multifloatsf)
        target_link_libraries(la_xisnan_mw PUBLIC
            $<BUILD_INTERFACE:multifloatsf>)
    endif()
    if(TARGET multifloats)
        target_link_libraries(la_xisnan_mw PUBLIC multifloats)
    endif()
endif()

# multifloats_mpi_f.f90: Fortran module exposing the C-side MPI
# datatype handles (MPI_FLOAT64X2 / MPI_DD_SUM / ...) via bind(c).
# MUMPS's migrated source `USE multifloats_mpi_f` requires the .mod;
# other libraries route MPI through C and don't need this target.
set(MF_MPI_DIR "{_mf_mpi_dir}"
    CACHE PATH "Directory containing multifloats_mpi_f.f90")
if(EXISTS "${{MF_MPI_DIR}}/multifloats_mpi_f.f90")
    add_library(multifloats_mpi_f STATIC
        "${{MF_MPI_DIR}}/multifloats_mpi_f.f90")
    set_target_properties(multifloats_mpi_f PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories(multifloats_mpi_f PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    if(MPI_Fortran_FOUND)
        target_link_libraries(multifloats_mpi_f PUBLIC MPI::MPI_Fortran)
    endif()
endif()
"""
            mf_deps = f"""
if(TARGET multifloats)
    target_link_libraries({precision_lib} PUBLIC multifloats)
endif()
# multifloatsf is the Fortran half of multifloats — provides the
# ``multifloats.mod`` module that the migrated source's ``use
# multifloats`` clauses resolve against. Wrap in $<BUILD_INTERFACE:>
# because multifloats owns its own install/export set; we just need
# the .mod path during this build.
if(TARGET multifloatsf)
    target_link_libraries({precision_lib} PUBLIC
        $<BUILD_INTERFACE:multifloatsf>)
endif()
if(TARGET la_constants_mw)
    target_link_libraries({precision_lib} PUBLIC la_constants_mw)
endif()
if(TARGET la_xisnan_mw)
    target_link_libraries({precision_lib} PUBLIC la_xisnan_mw)
endif()
if(TARGET multifloats_mpi_f)
    target_link_libraries({precision_lib} PUBLIC multifloats_mpi_f)
endif()
"""

        cmake = f"""\
cmake_minimum_required(VERSION 3.20)
project({precision_lib} Fortran C CXX)

# --- Compiler flags ---
set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -w")
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${{CMAKE_Fortran_FLAGS}} -std=legacy")
endif()

# Source-form line-length relief. The migrator lengthens tokens (e.g.
# MPI_DOUBLE_COMPLEX -> MPI_C_LONG_DOUBLE_COMPLEX), which can push
# fixed-form .F lines past column 72 and free-form .f90 lines past 132.
# Disable both limits, keyed on compiler family. Mirrors
# add_migrated_fortran_library in cmake/CMakeLists.txt (the canonical
# `migrator stage` path) so single-library `migrator build` builds stay
# consistent across gfortran / flang / Intel. CMAKE_Fortran_COMPILER_ID
# is the built-in analog of that file's FORTRAN_COMPILER_FAMILY.
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Flang")
    set(_fortran_line_length_flags
        $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(_fortran_line_length_flags
        $<$<COMPILE_LANGUAGE:Fortran>:-extend-source>)
else()
    set(_fortran_line_length_flags "")
endif()

# Enable Fortran preprocessing for .F90 files
set(CMAKE_Fortran_PREPROCESS ON)

# --- MPI: default to vendored Intel MPI headers ---
# See note in the C template: headers come from ``external/impi-headers/``
# unconditionally; the runtime comes from whichever MPI the user links
# against at final link time. MUMPS uses ``INCLUDE 'mpif.h'`` in 231
# source files and never ``USE mpi``, so F77 headers are enough.
if(NOT DEFINED IMPI_HEADERS)
    set(IMPI_HEADERS "{_impi_default}"
        CACHE PATH "Path to vendored Intel MPI headers")
endif()
include_directories(${{IMPI_HEADERS}})
find_package(MPI COMPONENTS Fortran QUIET)

# Detect extended-precision (KIND=10 / KIND=16) support.
# Shared probe sits in cmake/DetectExtendedPrecision.cmake; copied
# next to this CMakeLists.txt at generation time so the staging
# tree stays self-contained.
include(${{CMAKE_CURRENT_SOURCE_DIR}}/DetectExtendedPrecision.cmake)
{mf_link}
# --- Standard-precision sibling archive ---
# Built from upstream Fortran sources alongside the migrated archive.
# Carries the original S/D/C/Z entry points and the precision-
# independent helpers (LSAME, XERBLA, LA_XISNAN module, ...) that
# the migrated archive's bodies reference but don't ship themselves.
# The migrated archive PUBLIC-links this so downstreams resolve both
# symbol families through one link line. Modules also flow to the
# migrated build via the shared module directory (e.g. la_xisnan.mod
# from std-precision is what the la_xisnan_mw helper `use`s).
set(REF_SOURCES
    {ref_list}
)

if(REF_SOURCES)
    add_library({lib_name} STATIC ${{REF_SOURCES}})
    set_target_properties({lib_name} PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories({lib_name} PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    # Disable the line-length limit on the standard archive too — same
    # family-aware flags as the precision archive (set above) — so the
    # build stays consistent across the pair.
    target_compile_options({lib_name} PRIVATE ${{_fortran_line_length_flags}})
endif()

# --- Common (type-independent) library ---
set(COMMON_SOURCES
    {common_list}
)

# --- Precision-specific library ---
set(PRECISION_SOURCES
    {precision_list}
)

if(COMMON_SOURCES)
    add_library({common_lib} STATIC ${{COMMON_SOURCES}})
    set_target_properties({common_lib} PROPERTIES
        Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
    target_include_directories({common_lib} PUBLIC
        $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>)
    target_compile_options({common_lib} PRIVATE ${{_fortran_line_length_flags}})
endif()

add_library({precision_lib} STATIC ${{PRECISION_SOURCES}})
set_target_properties({precision_lib} PROPERTIES
    Fortran_MODULE_DIRECTORY ${{CMAKE_CURRENT_BINARY_DIR}}/mod)
target_include_directories({precision_lib} PUBLIC
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_BINARY_DIR}}/mod>
    $<BUILD_INTERFACE:${{CMAKE_CURRENT_SOURCE_DIR}}>)
target_compile_options({precision_lib} PRIVATE ${{_fortran_line_length_flags}})
if(TARGET {common_lib})
    target_link_libraries({precision_lib} PUBLIC {common_lib})
endif()
if(TARGET {lib_name})
    target_link_libraries({precision_lib} PUBLIC {lib_name})
endif()
if(MPI_Fortran_FOUND)
    target_link_libraries({precision_lib} PUBLIC MPI::MPI_Fortran)
    if(TARGET {common_lib})
        target_link_libraries({common_lib} PUBLIC MPI::MPI_Fortran)
    endif()
    if(TARGET {lib_name})
        target_link_libraries({lib_name} PUBLIC MPI::MPI_Fortran)
    endif()
endif()
{mf_deps}

# --- Install rules ---
install(TARGETS {precision_lib} ARCHIVE DESTINATION lib)
if(TARGET {common_lib})
    install(TARGETS {common_lib} ARCHIVE DESTINATION lib)
endif()
if(TARGET {lib_name})
    install(TARGETS {lib_name} ARCHIVE DESTINATION lib)
endif()
"""
    (output_dir / 'CMakeLists.txt').write_text(cmake)

    # Ship the shared extended-precision probe alongside the generated
    # CMakeLists.txt. The Fortran template ``include(...)``s it; the C
    # template doesn't need it but the file is cheap to copy and keeps
    # the staging tree self-contained.
    if language != 'c':
        _root = project_root or Path.cwd()
        probe = _root / 'cmake' / 'DetectExtendedPrecision.cmake'
        if probe.exists():
            shutil.copy2(probe, output_dir / probe.name)

