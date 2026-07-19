# FortranCompiler.cmake
#
# Detects the Fortran compiler's .mod file format version and provides
# functions for compiler-aware installation of Fortran modules and libraries.
#
# After include(FortranCompiler), the following variables are set:
#
#   FORTRAN_COMPILER_FAMILY   - Normalized compiler family (gfortran, intel, flang, nvhpc, nag, cray)
#   FORTRAN_COMPILER_VERSION  - Full compiler version string
#   FORTRAN_MOD_VERSION       - Internal .mod format version integer (or "unknown")
#   FORTRAN_MOD_COMPAT_TAG    - Module compat tag for .mod dirs (e.g. gfortran-mod15)
#   FORTRAN_COMPILER_TAG      - Compiler version tag for libraries (e.g. gfortran-14)
#
# Functions provided:
#
#   fortran_module_layout(<target>)
#     Sets up build-time and install-time module directories for the target.
#     Must be called before fortran_install_modules() or fortran_install_library().
#
#   fortran_install_modules(<target> [DESTINATION <base>])
#     Installs .mod/.smod files to <base>/fmod/<mod-tag>/.
#     Requires fortran_module_layout() to have been called on the target first.
#
#   fortran_install_library(<target> [NAMESPACE <ns>] [EXPORT <export-name>])
#     Installs the library with a compiler-tagged filename and generates a
#     Config.cmake for find_package() support.
#
#   _fc_detect_mpi_tag(<outvar>)
#     Probes mpi.h's vendor macros and returns the MPI implementation
#     flavor tag (e.g. intelmpi-2021.18), or "" when no MPI was found.
#     Called by the top-level CMakeLists to compute MPI_TAG; uses the
#     same canonical probe program the generated Configs re-run on the
#     consumer side, so the two tags cannot drift.

if(_FORTRAN_COMPILER_INCLUDED)
  return()
endif()
set(_FORTRAN_COMPILER_INCLUDED TRUE)

# ---------------------------------------------------------------------------
# Verify Fortran is enabled
# ---------------------------------------------------------------------------
get_property(_fc_languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if(NOT "Fortran" IN_LIST _fc_languages)
  message(FATAL_ERROR "FortranCompiler: Fortran language must be enabled before including this module.")
endif()
unset(_fc_languages)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---------------------------------------------------------------------------
# Detect compiler family
#
# Compiler IDs used by CMake:
#   GNU       - gfortran
#   Intel     - ifort (classic)
#   IntelLLVM - ifx
#   LLVMFlang - LLVM Flang (flang-new, the official LLVM Fortran compiler);
#               .mod files are valid Fortran source with !mod$ v1 header
#   Flang     - Classic Flang (PGI-derived, incompatible with LLVM Flang)
#   NVHPC     - NVIDIA nvfortran (PGI lineage, shares classic Flang .mod format)
#   NAG       - NAG Fortran
#   Cray      - Cray Fortran (CCE)
#
# The ID → family map and the per-family ABI version truncation are
# needed twice: once here at configure time (producer side, baked into
# the installed archive filenames) and once inside every generated
# <pkg>Config.cmake, where the consumer re-derives the same tag to pick
# the matching targets file. Both sides must stay in lockstep or
# find_package() stops matching the installed filenames — so the
# detection code exists exactly once, as text emitted by the two
# generators below, parameterized by output variable names. The
# producer evaluates it via cmake_language(EVAL CODE);
# fortran_install_library splices the identical text (with
# consumer-side names) into the generated Config.
# ---------------------------------------------------------------------------
function(_fc_family_detection_code outvar family_var)
  set(${outvar} "\
if(CMAKE_Fortran_COMPILER_ID STREQUAL \"GNU\")
  set(${family_var} \"gfortran\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Intel\")
  set(${family_var} \"intel\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"IntelLLVM\")
  set(${family_var} \"intel\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"LLVMFlang\")
  set(${family_var} \"flang\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Flang\")
  set(${family_var} \"flang-classic\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"NVHPC\")
  set(${family_var} \"nvhpc\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"NAG\")
  set(${family_var} \"nag\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Cray\")
  set(${family_var} \"cray\")
else()
  set(${family_var} \"\${CMAKE_Fortran_COMPILER_ID}\")
  string(TOLOWER \"\${${family_var}}\" ${family_var})
endif()
" PARENT_SCOPE)
endfunction()

# ABI-relevant version truncation per family, composing
# ``<family>-<abi-version>`` into ${tag_var} (and unsetting the scratch
# ${abiver_var}):
#   gfortran -> major only (ABI stable within a release series)
#   flang    -> major only (follows LLVM major versioning)
#   intel    -> major.minor (ABI can change at minor releases, e.g. 2021.10)
#   others   -> full version (conservative)
function(_fc_abi_tag_code outvar family_var version_var abiver_var tag_var)
  set(${outvar} "\
if(${family_var} STREQUAL \"gfortran\" OR ${family_var} STREQUAL \"flang\")
  string(REGEX MATCH \"^([0-9]+)\" ${abiver_var} \"\${${version_var}}\")
elseif(${family_var} STREQUAL \"intel\")
  string(REGEX MATCH \"^([0-9]+\\\\.[0-9]+)\" ${abiver_var} \"\${${version_var}}\")
else()
  set(${abiver_var} \"\${${version_var}}\")
endif()
set(${tag_var} \"\${${family_var}}-\${${abiver_var}}\")
unset(${abiver_var})
" PARENT_SCOPE)
endfunction()

set(FORTRAN_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")

_fc_family_detection_code(_fc_detect_code FORTRAN_COMPILER_FAMILY)
cmake_language(EVAL CODE "${_fc_detect_code}")
unset(_fc_detect_code)

# ---------------------------------------------------------------------------
# Determine .mod format version from compiler family + version
# ---------------------------------------------------------------------------
set(FORTRAN_MOD_VERSION "unknown")

if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran")
  # GCC major version -> MOD_VERSION
  # GCC < 4.4: unversioned .mod format (unsupported)
  string(REGEX MATCH "^([0-9]+)" _fc_gcc_major "${FORTRAN_COMPILER_VERSION}")
  if(_fc_gcc_major VERSION_GREATER_EQUAL 15)
    set(FORTRAN_MOD_VERSION "16")
  elseif(_fc_gcc_major VERSION_GREATER_EQUAL 8)
    set(FORTRAN_MOD_VERSION "15")
  elseif(_fc_gcc_major VERSION_GREATER_EQUAL 5)
    set(FORTRAN_MOD_VERSION "14")
  elseif(_fc_gcc_major EQUAL 4)
    string(REGEX MATCH "^4\\.([0-9]+)" _fc_gcc_4minor "${FORTRAN_COMPILER_VERSION}")
    set(_fc_minor "${CMAKE_MATCH_1}")
    if(_fc_minor EQUAL 9)
      set(FORTRAN_MOD_VERSION "12")
    elseif(_fc_minor EQUAL 8)
      set(FORTRAN_MOD_VERSION "10")
    elseif(_fc_minor EQUAL 7)
      set(FORTRAN_MOD_VERSION "9")
    elseif(_fc_minor EQUAL 6)
      set(FORTRAN_MOD_VERSION "6")
    elseif(_fc_minor EQUAL 5)
      set(FORTRAN_MOD_VERSION "4")
    elseif(_fc_minor EQUAL 4)
      set(FORTRAN_MOD_VERSION "0")
    endif()
    unset(_fc_minor)
    unset(_fc_gcc_4minor)
  endif()
  unset(_fc_gcc_major)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Classic ifort versioning
  # Version numbering jumped from 19.x to 2021.x (year-based oneAPI scheme).
  # All versions 18.x through 2021.x compare correctly with VERSION_GREATER_EQUAL
  # because 18 < 2020 < 2021 numerically.
  if(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2021.10")
    set(FORTRAN_MOD_VERSION "13")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "18.0")
    set(FORTRAN_MOD_VERSION "12")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "17.0")
    set(FORTRAN_MOD_VERSION "11")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "16.0")
    set(FORTRAN_MOD_VERSION "10")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  # ifx versioning (also uses year-based oneAPI scheme)
  if(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2023.2")
    set(FORTRAN_MOD_VERSION "13")
  elseif(FORTRAN_COMPILER_VERSION VERSION_GREATER_EQUAL "2021.0")
    set(FORTRAN_MOD_VERSION "12")
  endif()

elseif(FORTRAN_COMPILER_FAMILY STREQUAL "flang")
  # LLVM Flang: text-based .mod files with !mod$ v1 header
  set(FORTRAN_MOD_VERSION "1")

# flang-classic, nvhpc, nag, cray: no known .mod version mapping.
# They fall through to FORTRAN_MOD_VERSION = "unknown" and get tagged
# by full compiler version (conservative but safe).
endif()

# ---------------------------------------------------------------------------
# Build tags:
#   FORTRAN_MOD_COMPAT_TAG  - for .mod directories (by format compatibility)
#   FORTRAN_COMPILER_TAG    - for library files (by ABI-relevant version)
#
# Version truncation rules are documented on _fc_abi_tag_code above.
# ---------------------------------------------------------------------------
_fc_abi_tag_code(_fc_detect_code FORTRAN_COMPILER_FAMILY FORTRAN_COMPILER_VERSION _fc_abi_version FORTRAN_COMPILER_TAG)
cmake_language(EVAL CODE "${_fc_detect_code}")
unset(_fc_detect_code)

if(FORTRAN_MOD_VERSION STREQUAL "unknown")
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_TAG}")
else()
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_FAMILY}-mod${FORTRAN_MOD_VERSION}")
endif()

message(STATUS "FortranCompiler: compiler=${CMAKE_Fortran_COMPILER_ID} ${FORTRAN_COMPILER_VERSION}")
message(STATUS "FortranCompiler: family=${FORTRAN_COMPILER_FAMILY}, mod_version=${FORTRAN_MOD_VERSION}")
message(STATUS "FortranCompiler: mod_tag=${FORTRAN_MOD_COMPAT_TAG}, lib_tag=${FORTRAN_COMPILER_TAG}")

# ---------------------------------------------------------------------------
# MPI vendor probe (single source).
#
# One canonical C program identifies the MPI implementation flavor and
# version from mpi.h's vendor-specific macros:
#   Intel MPI : I_MPI_VERSION (string, e.g. "2021.18.0")
#   OpenMPI   : OMPI_MAJOR_VERSION / OMPI_MINOR_VERSION (ints)
#   MPICH     : MPICH_VERSION (string, e.g. "4.2.0")
# Its ``<flavor> <version>`` output parses via _FC_MPI_TAG_REGEX into a
# ``<flavor>-<major>.<minor>`` tag (``intelmpi-2021.18`` /
# ``openmpi-5.0`` / ``mpich-4.2``). The same program and regex serve
# the producer-side probe (_fc_detect_mpi_tag, called by the top-level
# CMakeLists to compute MPI_TAG) and the consumer-side probe re-escaped
# into every generated Config (_fc_consumer_mpi_detect_block below) —
# the two must derive identical tags or the installed filenames and the
# consumer's lookup diverge.
# ---------------------------------------------------------------------------
set(_FC_MPI_PROBE_SOURCE "\
#include <mpi.h>
#include <stdio.h>
int main(void){
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

# ``flavor major.minor[.patch]`` → capture flavor, major, minor.
set(_FC_MPI_TAG_REGEX "^([a-z]+) ([0-9]+)\\.([0-9]+)")

# _fc_detect_mpi_tag(<outvar>)
#
# Compile and run the canonical probe against the found MPI's headers
# and return ``<flavor>-<major>.<minor>``. Empty when no MPI was found,
# the probe fails, or the flavor is unidentified.
function(_fc_detect_mpi_tag outvar)
  set(_mpi_tag "")
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
    file(WRITE "${_mpi_probe_c}" "${_FC_MPI_PROBE_SOURCE}")
    try_run(_mpi_run_rc _mpi_compile_rc
        "${CMAKE_BINARY_DIR}/CMakeFiles/mpi_probe.dir"
        "${_mpi_probe_c}"
        CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${_mpi_inc}"
        RUN_OUTPUT_VARIABLE _mpi_id)
    if(_mpi_compile_rc AND "${_mpi_run_rc}" STREQUAL "0")
      string(STRIP "${_mpi_id}" _mpi_id)
      if(_mpi_id MATCHES "${_FC_MPI_TAG_REGEX}")
        set(_mpi_tag "${CMAKE_MATCH_1}-${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
      endif()
    endif()
  endif()
  set(${outvar} "${_mpi_tag}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# Helper: ensure INSTALL_INTERFACE paths are relative (required for
# relocatable packages). GNUInstallDirs can produce absolute paths on
# some platforms.
# ---------------------------------------------------------------------------
function(_fc_make_relative_installdir outvar path)
  if(IS_ABSOLUTE "${path}")
    file(RELATIVE_PATH _fc_relpath "${CMAKE_INSTALL_PREFIX}" "${path}")
    set(${outvar} "${_fc_relpath}" PARENT_SCOPE)
  else()
    set(${outvar} "${path}" PARENT_SCOPE)
  endif()
endfunction()

# ---------------------------------------------------------------------------
# fortran_module_layout(<target>)
#
# Configures the target's module output directory and include paths.
# Uses FORTRAN_MOD_COMPAT_TAG for module directories.
# Must be called before fortran_install_modules() or fortran_install_library().
# ---------------------------------------------------------------------------
function(fortran_module_layout target)
  set(_moddir "${PROJECT_BINARY_DIR}/fmod/${FORTRAN_MOD_COMPAT_TAG}")

  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY "${_moddir}"
  )

  # Ensure the install interface path is relative for relocatable packages
  _fc_make_relative_installdir(_fc_rel_libdir "${CMAKE_INSTALL_LIBDIR}")

  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${_moddir}>
    $<INSTALL_INTERFACE:${_fc_rel_libdir}/fmod/${FORTRAN_MOD_COMPAT_TAG}>
  )
endfunction()

# ---------------------------------------------------------------------------
# fortran_test_module_layout(<target>)
#
# Module layout for test-only targets. Modules go to a parallel
# fmod_tests/ tree instead of the shared fmod/ tree so that
# fortran_install_modules() — which installs the *entire* shared
# module directory — can never sweep test-harness modules into a
# release package. No INSTALL_INTERFACE path: test targets are never
# installed.
# ---------------------------------------------------------------------------
function(fortran_test_module_layout target)
  set(_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}")

  set_target_properties(${target} PROPERTIES
    Fortran_MODULE_DIRECTORY "${_moddir}"
  )

  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${_moddir}>
  )
endfunction()

# ---------------------------------------------------------------------------
# fortran_install_modules(<target> [DESTINATION <base>])
#
# Installs .mod and .smod files to <base>/fmod/<mod-tag>/
# Default DESTINATION is ${CMAKE_INSTALL_LIBDIR}.
# Requires fortran_module_layout() to have been called on the target first.
# ---------------------------------------------------------------------------
function(fortran_install_modules target)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "" "DESTINATION" "")
  if(ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "fortran_install_modules: unexpected arguments: ${ARG_UNPARSED_ARGUMENTS}")
  endif()

  # Validate that fortran_module_layout was called
  get_target_property(_fc_moddir ${target} Fortran_MODULE_DIRECTORY)
  if(NOT _fc_moddir)
    message(FATAL_ERROR
      "fortran_install_modules(${target}): Fortran_MODULE_DIRECTORY is not set. "
      "Call fortran_module_layout(${target}) first.")
  endif()

  if(NOT ARG_DESTINATION)
    set(ARG_DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  endif()

  set(_install_moddir "${ARG_DESTINATION}/fmod/${FORTRAN_MOD_COMPAT_TAG}")

  install(
    DIRECTORY "${_fc_moddir}/"
    DESTINATION "${_install_moddir}"
    FILES_MATCHING
      PATTERN "*.mod"
      PATTERN "*.smod"
  )
endfunction()

# ---------------------------------------------------------------------------
# Internal helpers for fortran_install_library
# ---------------------------------------------------------------------------

# TRUE iff <target> compiles Fortran sources (detected from the SOURCES
# property): only such archives depend on the Fortran compiler's name
# mangling, module format and runtime library. Pure-C archives have a
# stable platform ABI and carry no compiler tag.
function(_fc_target_has_fortran outvar target)
  get_target_property(_fc_sources ${target} SOURCES)
  set(_has_fortran FALSE)
  foreach(_fc_src IN LISTS _fc_sources)
    string(TOLOWER "${_fc_src}" _fc_src_lower)
    if(_fc_src_lower MATCHES "\\.(f|for|ftn|f77|f90|f95|f03|f08)$")
      set(_has_fortran TRUE)
      break()
    endif()
  endforeach()
  set(${outvar} "${_has_fortran}" PARENT_SCOPE)
endfunction()

# Assemble the install tag from the components the archive's ABI
# actually depends on — each appears only when genuinely needed.
# Compiler tag: only when <has_fortran> (name mangling, module format
# and runtime library differ per compiler). MPI flavor tag: only when
# <want_mpi> (objects are MPI-ABI-bound); prefer ``MPI_LIB_TAG`` (which
# a caller may override, e.g. to ``seq`` for the libmpiseq release) and
# fall back to the raw ``MPI_TAG``. When neither MPI variable is set,
# <want_mpi> is a no-op. Archives needing neither tag install untagged.
function(_fc_derive_install_tag out_tag out_is_mpi has_fortran want_mpi)
  set(_flavor_tag "${MPI_TAG}")
  if(DEFINED MPI_LIB_TAG)
    set(_flavor_tag "${MPI_LIB_TAG}")
  endif()

  set(_tag_parts "")
  if(has_fortran)
    list(APPEND _tag_parts "${FORTRAN_COMPILER_TAG}")
  endif()
  set(_is_mpi FALSE)
  if(want_mpi AND _flavor_tag)
    list(APPEND _tag_parts "${_flavor_tag}")
    set(_is_mpi TRUE)
  endif()
  list(JOIN _tag_parts "-" _full_tag)

  set(${out_tag} "${_full_tag}" PARENT_SCOPE)
  set(${out_is_mpi} "${_is_mpi}" PARENT_SCOPE)
endfunction()

# MPI language components (C;Fortran;CXX subset) whose MPI::MPI_*
# imported targets the exported Targets file will reference. Even a
# flavor-untagged package may reference them — static archives export
# their PRIVATE deps as $<LINK_ONLY:…>, and pblas/scalapack/xblas
# compile against mpi.h without their objects binding to the MPI ABI.
# Such a Config must find_dependency(MPI) or the imported-target
# references dangle at consumer time. Derive the needed components from
# the link interface rather than from the MPI flag. (MPI_C needs the
# delimiter guard so it doesn't match the MPI_CXX substring.)
function(_fc_mpi_dep_components outvar target)
  get_target_property(_fc_iface_links ${target} INTERFACE_LINK_LIBRARIES)
  set(_components "")
  if(_fc_iface_links MATCHES "MPI::MPI_C(>|;|$)")
    list(APPEND _components C)
  endif()
  if(_fc_iface_links MATCHES "MPI::MPI_Fortran")
    list(APPEND _components Fortran)
  endif()
  if(_fc_iface_links MATCHES "MPI::MPI_CXX")
    list(APPEND _components CXX)
  endif()
  set(${outvar} "${_components}" PARENT_SCOPE)
endfunction()

# Consumer-side compiler detection block spliced into the Config of
# every archive with Fortran objects. Pure-C packages skip it entirely,
# so a consumer project with no Fortran language enabled can still
# find_package() them. The detection code itself comes from the same
# generators the producer-side detection evaluated above — only the
# variable names differ.
function(_fc_consumer_fc_detect_block outvar)
  _fc_family_detection_code(_family_code _FC_consumer_family)
  _fc_abi_tag_code(_abi_code _FC_consumer_family CMAKE_Fortran_COMPILER_VERSION _FC_abi_version _FC_consumer_tag)
  set(${outvar} "\
# --- Derive consumer's compiler family and ABI version tag ---
set(_FC_consumer_family \"\")

${_family_code}
${_abi_code}" PARENT_SCOPE)
endfunction()

# Consumer-side MPI flavor detection block spliced into the Config of
# MPI-ABI-bound (flavor-tagged) archives. Re-runs the canonical
# _FC_MPI_PROBE_SOURCE program against the consumer's mpi.h.
function(_fc_consumer_mpi_detect_block outvar config_name)
  # Escape the canonical probe program and tag regex to survive one
  # level of CMake string parsing inside the emitted Config.
  string(REPLACE "\\" "\\\\" _probe "${_FC_MPI_PROBE_SOURCE}")
  string(REPLACE "\"" "\\\"" _probe "${_probe}")
  string(REPLACE "#" "\\#" _probe "${_probe}")
  string(REPLACE "\n" "\\n" _probe "${_probe}")
  string(REPLACE "\\" "\\\\" _regex "${_FC_MPI_TAG_REGEX}")
  set(${outvar} "\
# --- Derive consumer's MPI flavor + major.minor tag ---
find_package(MPI QUIET COMPONENTS C)
set(_FC_consumer_mpi_tag \"\")
# A consumer may FORCE the MPI flavor tag via -DEPLINALG_MPI_TAG=<tag>,
# bypassing the mpi.h vendor probe below. This is required for a
# libmpiseq/seq release: such a consumer still finds a real MPI for the
# mpi.h *headers* (mmsolve.c #includes mpi.h), so the probe would detect
# that vendor (e.g. intelmpi-2021.18) — but the installed archives are
# tagged `seq`. Setting EPLINALG_MPI_TAG=seq resolves the seq targets file.
if(DEFINED EPLINALG_MPI_TAG)
  set(_FC_consumer_mpi_tag \"\${EPLINALG_MPI_TAG}\")
elseif(MPI_C_FOUND)
  set(_FC_mpi_inc \"\${MPI_C_HEADER_DIR}\")
  if(TARGET MPI::MPI_C)
    get_target_property(_FC_mpi_iface MPI::MPI_C INTERFACE_INCLUDE_DIRECTORIES)
    if(_FC_mpi_iface)
      list(APPEND _FC_mpi_inc \${_FC_mpi_iface})
    endif()
  endif()
  list(REMOVE_DUPLICATES _FC_mpi_inc)
  set(_FC_mpi_probe \"\${CMAKE_CURRENT_BINARY_DIR}/_FC_mpi_probe.c\")
  file(WRITE \"\${_FC_mpi_probe}\" \"${_probe}\")
  try_run(_FC_run_rc _FC_compile_rc
    \"\${CMAKE_CURRENT_BINARY_DIR}/_FC_mpi_probe.dir\"
    \"\${_FC_mpi_probe}\"
    CMAKE_FLAGS \"-DINCLUDE_DIRECTORIES=\${_FC_mpi_inc}\"
    RUN_OUTPUT_VARIABLE _FC_mpi_id)
  if(_FC_compile_rc AND \"\${_FC_run_rc}\" STREQUAL \"0\")
    string(STRIP \"\${_FC_mpi_id}\" _FC_mpi_id)
    if(_FC_mpi_id MATCHES \"${_regex}\")
      set(_FC_consumer_mpi_tag \"\${CMAKE_MATCH_1}-\${CMAKE_MATCH_2}.\${CMAKE_MATCH_3}\")
    endif()
  endif()
  unset(_FC_mpi_inc)
  unset(_FC_mpi_iface)
  unset(_FC_mpi_probe)
  unset(_FC_mpi_id)
  unset(_FC_run_rc)
  unset(_FC_compile_rc)
endif()
if(NOT _FC_consumer_mpi_tag)
  set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    \"${config_name}: MPI-dependent library but no MPI flavor detected on the consumer side. find_package(MPI) must succeed and mpi.h must be locatable.\")
  return()
endif()
" PARENT_SCOPE)
endfunction()

# Content of the per-tag transparent-deps file installed next to the
# per-tag targets file (see the flavor comment at the call site).
function(_fc_tag_deps_content outvar deps_file config_name full_tag mpi_dep_block deps_block)
  set(${outvar} "\
# ${deps_file}
# Auto-generated by FortranCompiler.cmake
#
# Transparent dependencies of the '${full_tag}' flavor of package
# ${config_name}. Kept per-tag (not in the shared Config) because one
# prefix can hold several flavors of this package with different
# dependency sets. Included by ${config_name}Config.cmake after it
# resolves the consumer's tag; CMakeFindDependencyMacro is already in
# scope there.
${mpi_dep_block}${deps_block}" PARENT_SCOPE)
endfunction()

# Content of the Config for a tagged package: re-derives the tag on the
# consumer side (via the detection blocks passed in) and includes the
# matching per-tag deps + targets files.
function(_fc_tagged_config_content outvar config_name export_name consumer_tag_expr fc_detect_block mpi_detect_block cleanup_mpi)
  set(${outvar} "\
# ${config_name}Config.cmake
# Auto-generated by FortranCompiler.cmake
#
# Re-derives, on the consumer side, the ABI tag components baked into
# this package's archive filename (compiler tag → Fortran ABI, MPI
# flavor tag → MPI ABI) and includes the matching targets file.
# Module directories are tagged by .mod format version separately.
#
# This file is identical across all flavors of the package and may be
# overwritten by any of them; everything flavor-specific (the archive
# targets AND their transparent find_dependency() calls) lives in the
# per-tag targets/deps files included below.

cmake_minimum_required(VERSION 3.12)

# find_dependency() is used by the per-tag deps file included below (the
# MPI dependency line when the Targets file references MPI targets, and
# the transparent DEPENDS block for libraries that link a factored-out
# shared package). No precision-sibling find_dependency() calls are
# emitted — consumers list every per-precision package they need on
# their own.
include(CMakeFindDependencyMacro)
${fc_detect_block}
${mpi_detect_block}
# Look for exact tag match
set(_FC_targets_file \"\${CMAKE_CURRENT_LIST_DIR}/${export_name}-${consumer_tag_expr}.cmake\")

if(NOT EXISTS \"\${_FC_targets_file}\")
  # No exact match — list available builds and fail
  file(GLOB _FC_available \"\${CMAKE_CURRENT_LIST_DIR}/${export_name}-*.cmake\")
  set(_FC_available_names \"\")
  foreach(_FC_f IN LISTS _FC_available)
    get_filename_component(_FC_fname \"\${_FC_f}\" NAME)
    list(APPEND _FC_available_names \"\${_FC_fname}\")
  endforeach()
  list(JOIN _FC_available_names \", \" _FC_available_list)
  set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    \"${config_name}: no pre-built library found for tag '${consumer_tag_expr}'. Available: [\${_FC_available_list}]\")
  unset(_FC_consumer_family)
  unset(_FC_consumer_tag)
  unset(_FC_targets_file)
  unset(_FC_available)
  unset(_FC_available_names)
  unset(_FC_available_list)
  ${cleanup_mpi}
  return()
endif()

# Transparent deps of THIS flavor (other flavors in the same prefix may
# have different ones), then the flavor's targets. The deps file runs
# nested find_dependency() calls whose Configs come from this same
# template and (re)use the same _FC_* variables in this scope — stash
# the targets path in a package-unique variable across that include.
set(_FC_targets_file_${export_name} \"\${_FC_targets_file}\")
include(\"\${CMAKE_CURRENT_LIST_DIR}/${export_name}-deps-${consumer_tag_expr}.cmake\")
include(\"\${_FC_targets_file_${export_name}}\")
unset(_FC_targets_file_${export_name})

unset(_FC_consumer_family)
unset(_FC_consumer_tag)
unset(_FC_targets_file)
${cleanup_mpi}
" PARENT_SCOPE)
endfunction()

# Content of the Config for an untagged package: a single fixed targets
# file, no consumer-side detection at all.
function(_fc_untagged_config_content outvar config_name export_name mpi_dep_block deps_block)
  set(${outvar} "\
# ${config_name}Config.cmake
# Auto-generated by FortranCompiler.cmake
#
# This package's archive contains no Fortran objects and is not bound
# to an MPI ABI, so it is compatible across Fortran compilers and MPI
# flavors — a single untagged targets file suffices and no
# consumer-side tag detection is performed.

cmake_minimum_required(VERSION 3.12)

include(CMakeFindDependencyMacro)
${mpi_dep_block}${deps_block}
include(\"\${CMAKE_CURRENT_LIST_DIR}/${export_name}.cmake\")
" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# fortran_install_library(<target>
#     [MPI]
#     [NAMESPACE <ns>]
#     [EXPORT <export-name>]
#     [DESTINATION <lib-dir>])
#
# Installs the library with a filename tagged by exactly the ABI axes the
# archive depends on, while the export/config system uses the mod compat
# tag to find modules:
#
#   - compiler tag (``gfortran-13``, …): included only when the target
#     compiles Fortran sources (auto-detected from the SOURCES property).
#     Pure-C archives have a stable platform ABI and carry no compiler tag.
#   - MPI flavor tag (``openmpi-4.1`` / ``mpich-4.2`` / ``seq``, …):
#     included only when ``MPI`` is passed, i.e. when the objects
#     themselves are MPI-ABI-bound (reference MPI symbols or bake mpi.h
#     constants). The flavor comes from ``MPI_LIB_TAG`` when the caller
#     defines it (letting a libmpiseq release override it to ``seq``),
#     otherwise from the raw ``MPI_TAG`` the top-level CMakeLists detects
#     from mpi.h's vendor macros. When neither is set, ``MPI`` is a no-op.
#
# An archive needing neither tag installs untagged (``libfoo.a``), and its
# Config performs no consumer-side detection at all. Independently of the
# flavor tag, a Config emits ``find_dependency(MPI COMPONENTS …)`` whenever
# the target's link interface references ``MPI::MPI_*`` imported targets
# (static archives export PRIVATE deps as ``$<LINK_ONLY:…>``), so e.g. an
# MPI-ABI-clean archive that merely compiles against mpi.h still resolves.
#
# Note: This function generates a <ProjectName>Config.cmake. If your project
# has multiple Fortran library targets, add them all to a single EXPORT set
# by passing the same EXPORT name. The Config.cmake and export file are
# generated only once per export set.
# ---------------------------------------------------------------------------
function(fortran_install_library target)
  cmake_parse_arguments(PARSE_ARGV 1 ARG "MPI;KEEP_OUTPUT_NAME" "NAMESPACE;EXPORT;DESTINATION;OUTPUT_BASE" "DEPENDS")
  if(ARG_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "fortran_install_library: unexpected arguments: ${ARG_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT ARG_NAMESPACE)
    set(ARG_NAMESPACE "${PROJECT_NAME}::")
  endif()
  if(NOT ARG_EXPORT)
    set(ARG_EXPORT "${PROJECT_NAME}Targets")
  endif()
  if(NOT ARG_DESTINATION)
    set(ARG_DESTINATION "${CMAKE_INSTALL_LIBDIR}")
  endif()

  _fc_target_has_fortran(_has_fortran ${target})
  _fc_derive_install_tag(_full_tag _is_mpi_lib "${_has_fortran}" "${ARG_MPI}")
  _fc_mpi_dep_components(_mpi_dep_components ${target})

  # Derive config name from the export set name (strip trailing "Targets").
  # This allows each library to get its own Config.cmake when given a
  # unique EXPORT name (e.g. EXPORT qblasTargets → qblasConfig.cmake).
  string(REGEX REPLACE "Targets$" "" _config_name "${ARG_EXPORT}")
  set(_cmake_install_dir "${ARG_DESTINATION}/cmake/${_config_name}")

  # Tag the library output filename. The target name IS the archive
  # base: migrated precision targets are pair-prefixed (eylapack →
  # libeylapack-<tag>.a), so name, filename and package agree. An
  # untagged target keeps its default output name (libeypblas.a).
  # OUTPUT_BASE substitutes a different archive base while the target
  # keeps its generic name (the ep_-privatized ``*_common`` archives
  # install as libep<lib>_common-<tag>.a but stay ``<lib>_common`` in
  # CMake). KEEP_OUTPUT_NAME preserves an OUTPUT_NAME already baked at
  # target creation (vendored archives like libptscotch_mumps-<tag>.a
  # carry a _mumps suffix plus their own flavor tag) — only the cmake
  # package machinery applies, not the rename.
  set(_output_base "${target}")
  if(ARG_OUTPUT_BASE)
    set(_output_base "${ARG_OUTPUT_BASE}")
  endif()
  if(_full_tag)
    set(_targets_file "${ARG_EXPORT}-${_full_tag}.cmake")
    if(NOT ARG_KEEP_OUTPUT_NAME)
      set_target_properties(${target} PROPERTIES
        OUTPUT_NAME "${_output_base}-${_full_tag}"
      )
    endif()
  else()
    set(_targets_file "${ARG_EXPORT}.cmake")
    if(ARG_OUTPUT_BASE AND NOT ARG_KEEP_OUTPUT_NAME)
      set_target_properties(${target} PROPERTIES
        OUTPUT_NAME "${_output_base}"
      )
    endif()
  endif()

  # Add target to the export set
  install(TARGETS ${target}
    EXPORT ${ARG_EXPORT}
    ARCHIVE DESTINATION "${ARG_DESTINATION}"
    LIBRARY DESTINATION "${ARG_DESTINATION}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )

  # Only generate the export file and Config.cmake once per export set.
  # Multiple targets can share the same EXPORT name; CMake accumulates them
  # from all install(TARGETS ... EXPORT <name>) calls. But install(EXPORT)
  # and Config.cmake generation must happen exactly once.
  get_property(_fc_config_generated GLOBAL PROPERTY _FC_CONFIG_GENERATED_${ARG_EXPORT})
  if(_fc_config_generated)
    return()
  endif()
  set_property(GLOBAL PROPERTY _FC_CONFIG_GENERATED_${ARG_EXPORT} TRUE)

  # Write the export set (includes all targets added to this EXPORT name)
  install(EXPORT ${ARG_EXPORT}
    FILE "${_targets_file}"
    NAMESPACE "${ARG_NAMESPACE}"
    DESTINATION "${_cmake_install_dir}"
  )

  # No cross-package sibling find_dependency() calls are emitted
  # into the generated Config.cmake — consumers list every archive
  # they need on their own link line; symbol resolution happens at
  # final link. Rationale: precision siblings (qxblas vs eyblas, …) are
  # independently useful, and auto-loading them would force every
  # consumer to install all of them. Only transparent deps a package's
  # own Targets file references (DEPENDS below) are emitted.

  # Generate Config.cmake that finds the right targets file. Strategy:
  # re-derive, on the consumer side, exactly the tag components this
  # install baked into the filename — compiler tag from the consumer's
  # Fortran compiler, MPI flavor tag from an mpi.h vendor probe — and
  # look for an exact match, failing with an informative error listing
  # available builds. A package with no tag components includes a
  # fixed targets file with no detection at all.

  # find_dependency(MPI) whenever the Targets file references
  # MPI::MPI_* — independent of the flavor tag (see
  # _fc_mpi_dep_components above).
  if(_mpi_dep_components)
    list(JOIN _mpi_dep_components " " _mpi_dep_components_str)
    set(_mpi_dep_block "\
# The Targets file references MPI::MPI_* imported targets (static
# archives export PRIVATE deps as LINK_ONLY entries); bring them into
# scope so those references resolve.
find_dependency(MPI COMPONENTS ${_mpi_dep_components_str})
")
  else()
    set(_mpi_dep_block "")
  endif()

  if(_is_mpi_lib)
    _fc_consumer_mpi_detect_block(_mpi_detect_block "${_config_name}")
    set(_cleanup_mpi "unset(_FC_consumer_mpi_tag)")
  else()
    set(_mpi_detect_block "")
    set(_cleanup_mpi "")
  endif()

  # Consumer-side tag expression mirrors the producer-side _tag_parts
  # composition: compiler part iff the archive has Fortran objects,
  # flavor part iff it is MPI-ABI-bound.
  set(_consumer_tag_parts "")
  if(_has_fortran)
    list(APPEND _consumer_tag_parts "\${_FC_consumer_tag}")
  endif()
  if(_is_mpi_lib)
    list(APPEND _consumer_tag_parts "\${_FC_consumer_mpi_tag}")
  endif()
  list(JOIN _consumer_tag_parts "-" _consumer_tag_expr)

  # Build a block of find_dependency() calls for the DEPENDS list. These
  # are transparent dependencies (not user-facing): factored-out shared
  # packages that the precision archive PUBLIC-links and that are each
  # installed as their own Config. Examples: the standard-precision
  # archive `eplinalg::blas` (package `eplinalgStdBlas`), and the
  # first-party MPI bridges (multifloats_mpi, quad_mpi). The
  # per-precision Config auto-loads them so consumers only need to
  # find_package(qxblas) / find_package(mwblas).
  set(_deps_block "")
  foreach(_dep IN LISTS ARG_DEPENDS)
    string(APPEND _deps_block "find_dependency(${_dep})\n")
  endforeach()

  if(_has_fortran)
    _fc_consumer_fc_detect_block(_fc_detect_block)
  else()
    set(_fc_detect_block "")
  endif()

  if(_full_tag)
    # The dependency set can differ BETWEEN FLAVORS of the same package
    # sharing one install prefix (e.g. a real-MPI mumps_common transparently
    # depends on PT-Scotch; a seq (libmpiseq) one does not, and its Targets
    # file references eplinalg::mpiseq instead of MPI::MPI_*). The Config
    # file below is flavor-independent and gets overwritten by whichever
    # flavor installs last — so the find_dependency() calls must NOT live in
    # it. Split them into a per-tag deps file installed next to the per-tag
    # targets file; the Config includes the one matching the consumer's tag.
    set(_deps_file "${ARG_EXPORT}-deps-${_full_tag}.cmake")
    _fc_tag_deps_content(_deps_content
      "${_deps_file}" "${_config_name}" "${_full_tag}"
      "${_mpi_dep_block}" "${_deps_block}")
    file(GENERATE
      OUTPUT "${PROJECT_BINARY_DIR}/cmake/${_deps_file}"
      CONTENT "${_deps_content}"
    )
    install(
      FILES "${PROJECT_BINARY_DIR}/cmake/${_deps_file}"
      DESTINATION "${_cmake_install_dir}"
    )

    _fc_tagged_config_content(_config_content
      "${_config_name}" "${ARG_EXPORT}" "${_consumer_tag_expr}"
      "${_fc_detect_block}" "${_mpi_detect_block}" "${_cleanup_mpi}")
  else()
    _fc_untagged_config_content(_config_content
      "${_config_name}" "${ARG_EXPORT}"
      "${_mpi_dep_block}" "${_deps_block}")
  endif()

  file(GENERATE
    OUTPUT "${PROJECT_BINARY_DIR}/cmake/${_config_name}Config.cmake"
    CONTENT "${_config_content}"
  )

  install(
    FILES "${PROJECT_BINARY_DIR}/cmake/${_config_name}Config.cmake"
    DESTINATION "${_cmake_install_dir}"
  )
endfunction()
