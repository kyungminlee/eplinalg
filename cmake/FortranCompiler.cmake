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
#   LLVMFlang - LLVM Flang (flang-new, the official LLVM Fortran compiler)
#   Flang     - Classic Flang (PGI-derived, incompatible with LLVM Flang)
#   NVHPC     - NVIDIA nvfortran (PGI lineage, shares classic Flang .mod format)
#   NAG       - NAG Fortran
#   Cray      - Cray Fortran (CCE)
# ---------------------------------------------------------------------------
set(FORTRAN_COMPILER_VERSION "${CMAKE_Fortran_COMPILER_VERSION}")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(FORTRAN_COMPILER_FAMILY "gfortran")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(FORTRAN_COMPILER_FAMILY "intel")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
  set(FORTRAN_COMPILER_FAMILY "intel")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "LLVMFlang")
  # LLVM Flang (flang-new): .mod files are valid Fortran source with !mod$ v1 header
  set(FORTRAN_COMPILER_FAMILY "flang")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Flang")
  # Classic Flang (PGI-derived): completely different .mod format from LLVM Flang
  set(FORTRAN_COMPILER_FAMILY "flang-classic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  # NVIDIA nvfortran: shares PGI-lineage .mod format with classic Flang
  set(FORTRAN_COMPILER_FAMILY "nvhpc")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
  set(FORTRAN_COMPILER_FAMILY "nag")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  set(FORTRAN_COMPILER_FAMILY "cray")
else()
  set(FORTRAN_COMPILER_FAMILY "${CMAKE_Fortran_COMPILER_ID}")
  string(TOLOWER "${FORTRAN_COMPILER_FAMILY}" FORTRAN_COMPILER_FAMILY)
endif()

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
# Version truncation per family:
#   gfortran -> major only (ABI stable within a release series)
#   flang    -> major only (follows LLVM major versioning)
#   intel    -> major.minor (ABI can change at minor releases, e.g. 2021.10)
#   others   -> full version (conservative)
# ---------------------------------------------------------------------------
if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran" OR FORTRAN_COMPILER_FAMILY STREQUAL "flang")
  string(REGEX MATCH "^([0-9]+)" _fc_abi_version "${FORTRAN_COMPILER_VERSION}")
elseif(FORTRAN_COMPILER_FAMILY STREQUAL "intel")
  string(REGEX MATCH "^([0-9]+\\.[0-9]+)" _fc_abi_version "${FORTRAN_COMPILER_VERSION}")
else()
  set(_fc_abi_version "${FORTRAN_COMPILER_VERSION}")
endif()
set(FORTRAN_COMPILER_TAG "${FORTRAN_COMPILER_FAMILY}-${_fc_abi_version}")
unset(_fc_abi_version)

if(FORTRAN_MOD_VERSION STREQUAL "unknown")
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_TAG}")
else()
  set(FORTRAN_MOD_COMPAT_TAG "${FORTRAN_COMPILER_FAMILY}-mod${FORTRAN_MOD_VERSION}")
endif()

message(STATUS "FortranCompiler: compiler=${CMAKE_Fortran_COMPILER_ID} ${FORTRAN_COMPILER_VERSION}")
message(STATUS "FortranCompiler: family=${FORTRAN_COMPILER_FAMILY}, mod_version=${FORTRAN_MOD_VERSION}")
message(STATUS "FortranCompiler: mod_tag=${FORTRAN_MOD_COMPAT_TAG}, lib_tag=${FORTRAN_COMPILER_TAG}")

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
  cmake_parse_arguments(PARSE_ARGV 1 ARG "MPI;KEEP_OUTPUT_NAME" "NAMESPACE;EXPORT;DESTINATION" "DEPENDS")
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

  # Assemble the install tag from the components the archive's ABI
  # actually depends on — each appears only when genuinely needed.
  # Compiler tag: only for targets that compile Fortran (name mangling,
  # module format and runtime library differ per compiler); detected
  # from the SOURCES property. MPI flavor tag: only when the caller
  # passes ``MPI`` (objects are MPI-ABI-bound); prefer ``MPI_LIB_TAG``
  # (which a caller may override, e.g. to ``seq`` for the libmpiseq
  # release) and fall back to the raw ``MPI_TAG``. Archives needing
  # neither install untagged.
  get_target_property(_fc_sources ${target} SOURCES)
  set(_has_fortran FALSE)
  foreach(_fc_src IN LISTS _fc_sources)
    string(TOLOWER "${_fc_src}" _fc_src_lower)
    if(_fc_src_lower MATCHES "\\.(f|for|ftn|f77|f90|f95|f03|f08)$")
      set(_has_fortran TRUE)
      break()
    endif()
  endforeach()

  set(_flavor_tag "${MPI_TAG}")
  if(DEFINED MPI_LIB_TAG)
    set(_flavor_tag "${MPI_LIB_TAG}")
  endif()

  set(_tag_parts "")
  if(_has_fortran)
    list(APPEND _tag_parts "${FORTRAN_COMPILER_TAG}")
  endif()
  set(_is_mpi_lib FALSE)
  if(ARG_MPI AND _flavor_tag)
    list(APPEND _tag_parts "${_flavor_tag}")
    set(_is_mpi_lib TRUE)
  endif()
  list(JOIN _tag_parts "-" _full_tag)

  # Even a flavor-untagged package may reference MPI::MPI_* imported
  # targets in its exported Targets file — static archives export their
  # PRIVATE deps as $<LINK_ONLY:…>, and pblas/scalapack/xblas compile
  # against mpi.h without their objects binding to the MPI ABI. Such a
  # Config must find_dependency(MPI) or the imported-target references
  # dangle at consumer time. Derive the needed components from the link
  # interface rather than from the MPI flag. (MPI_C needs the delimiter
  # guard so it doesn't match the MPI_CXX substring.)
  get_target_property(_fc_iface_links ${target} INTERFACE_LINK_LIBRARIES)
  set(_mpi_dep_components "")
  if(_fc_iface_links MATCHES "MPI::MPI_C(>|;|$)")
    list(APPEND _mpi_dep_components C)
  endif()
  if(_fc_iface_links MATCHES "MPI::MPI_Fortran")
    list(APPEND _mpi_dep_components Fortran)
  endif()
  if(_fc_iface_links MATCHES "MPI::MPI_CXX")
    list(APPEND _mpi_dep_components CXX)
  endif()

  # Derive config name from the export set name (strip trailing "Targets").
  # This allows each library to get its own Config.cmake when given a
  # unique EXPORT name (e.g. EXPORT qblasTargets → qblasConfig.cmake).
  string(REGEX REPLACE "Targets$" "" _config_name "${ARG_EXPORT}")
  set(_cmake_install_dir "${ARG_DESTINATION}/cmake/${_config_name}")

  # Tag the library output filename. The target name IS the archive
  # base: migrated precision targets are pair-prefixed (eylapack →
  # libeylapack-<tag>.a), so name, filename and package agree. An
  # untagged target keeps its default output name (libeypblas.a).
  # KEEP_OUTPUT_NAME preserves an OUTPUT_NAME already baked at target
  # creation (vendored archives like libptscotch_mumps-<tag>.a carry a
  # _mumps suffix plus their own flavor tag) — only the cmake package
  # machinery applies, not the rename.
  if(_full_tag)
    set(_targets_file "${ARG_EXPORT}-${_full_tag}.cmake")
    if(NOT ARG_KEEP_OUTPUT_NAME)
      set_target_properties(${target} PROPERTIES
        OUTPUT_NAME "${target}-${_full_tag}"
      )
    endif()
  else()
    set(_targets_file "${ARG_EXPORT}.cmake")
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
  # MPI::MPI_* — independent of the flavor tag (see the derivation of
  # _mpi_dep_components above).
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
    set(_mpi_detect_block "\
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
  file(WRITE \"\${_FC_mpi_probe}\" \"\\#include <mpi.h>\\n\\#include <stdio.h>\\nint main(void){\\n\\#if defined(I_MPI_VERSION)\\n  printf(\\\"intelmpi %s\\\\n\\\", I_MPI_VERSION);\\n\\#elif defined(OMPI_MAJOR_VERSION) && defined(OMPI_MINOR_VERSION)\\n  printf(\\\"openmpi %d.%d\\\\n\\\", OMPI_MAJOR_VERSION, OMPI_MINOR_VERSION);\\n\\#elif defined(MPICH_VERSION)\\n  printf(\\\"mpich %s\\\\n\\\", MPICH_VERSION);\\n\\#else\\n  printf(\\\"unknown ?\\\\n\\\");\\n\\#endif\\n  return 0;\\n}\\n\")
  try_run(_FC_run_rc _FC_compile_rc
    \"\${CMAKE_CURRENT_BINARY_DIR}/_FC_mpi_probe.dir\"
    \"\${_FC_mpi_probe}\"
    CMAKE_FLAGS \"-DINCLUDE_DIRECTORIES=\${_FC_mpi_inc}\"
    RUN_OUTPUT_VARIABLE _FC_mpi_id)
  if(_FC_compile_rc AND \"\${_FC_run_rc}\" STREQUAL \"0\")
    string(STRIP \"\${_FC_mpi_id}\" _FC_mpi_id)
    if(_FC_mpi_id MATCHES \"^([a-z]+) ([0-9]+)\\\\.([0-9]+)\")
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
    \"${_config_name}: MPI-dependent library but no MPI flavor detected on the consumer side. find_package(MPI) must succeed and mpi.h must be locatable.\")
  return()
endif()
")
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

  # Compiler-detection block: emitted only for archives with Fortran
  # objects. Pure-C packages skip it entirely, so a consumer project
  # with no Fortran language enabled can still find_package() them.
  if(_has_fortran)
    set(_fc_detect_block "\
# --- Derive consumer's compiler family and ABI version tag ---
set(_FC_consumer_family \"\")

if(CMAKE_Fortran_COMPILER_ID STREQUAL \"GNU\")
  set(_FC_consumer_family \"gfortran\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Intel\")
  set(_FC_consumer_family \"intel\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"IntelLLVM\")
  set(_FC_consumer_family \"intel\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"LLVMFlang\")
  set(_FC_consumer_family \"flang\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Flang\")
  set(_FC_consumer_family \"flang-classic\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"NVHPC\")
  set(_FC_consumer_family \"nvhpc\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"NAG\")
  set(_FC_consumer_family \"nag\")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL \"Cray\")
  set(_FC_consumer_family \"cray\")
else()
  set(_FC_consumer_family \"\${CMAKE_Fortran_COMPILER_ID}\")
  string(TOLOWER \"\${_FC_consumer_family}\" _FC_consumer_family)
endif()

if(_FC_consumer_family STREQUAL \"gfortran\" OR _FC_consumer_family STREQUAL \"flang\")
  string(REGEX MATCH \"^([0-9]+)\" _FC_abi_version \"\${CMAKE_Fortran_COMPILER_VERSION}\")
elseif(_FC_consumer_family STREQUAL \"intel\")
  string(REGEX MATCH \"^([0-9]+\\\\.[0-9]+)\" _FC_abi_version \"\${CMAKE_Fortran_COMPILER_VERSION}\")
else()
  set(_FC_abi_version \"\${CMAKE_Fortran_COMPILER_VERSION}\")
endif()
set(_FC_consumer_tag \"\${_FC_consumer_family}-\${_FC_abi_version}\")
unset(_FC_abi_version)
")
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
    file(GENERATE
      OUTPUT "${PROJECT_BINARY_DIR}/cmake/${_deps_file}"
      CONTENT "\
# ${_deps_file}
# Auto-generated by FortranCompiler.cmake
#
# Transparent dependencies of the '${_full_tag}' flavor of package
# ${_config_name}. Kept per-tag (not in the shared Config) because one
# prefix can hold several flavors of this package with different
# dependency sets. Included by ${_config_name}Config.cmake after it
# resolves the consumer's tag; CMakeFindDependencyMacro is already in
# scope there.
${_mpi_dep_block}${_deps_block}"
    )
    install(
      FILES "${PROJECT_BINARY_DIR}/cmake/${_deps_file}"
      DESTINATION "${_cmake_install_dir}"
    )

    set(_config_content "\
# ${_config_name}Config.cmake
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
${_fc_detect_block}
${_mpi_detect_block}
# Look for exact tag match
set(_FC_targets_file \"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}-${_consumer_tag_expr}.cmake\")

if(NOT EXISTS \"\${_FC_targets_file}\")
  # No exact match — list available builds and fail
  file(GLOB _FC_available \"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}-*.cmake\")
  set(_FC_available_names \"\")
  foreach(_FC_f IN LISTS _FC_available)
    get_filename_component(_FC_fname \"\${_FC_f}\" NAME)
    list(APPEND _FC_available_names \"\${_FC_fname}\")
  endforeach()
  list(JOIN _FC_available_names \", \" _FC_available_list)
  set(\${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
  set(\${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE
    \"${_config_name}: no pre-built library found for tag '${_consumer_tag_expr}'. Available: [\${_FC_available_list}]\")
  unset(_FC_consumer_family)
  unset(_FC_consumer_tag)
  unset(_FC_targets_file)
  unset(_FC_available)
  unset(_FC_available_names)
  unset(_FC_available_list)
  ${_cleanup_mpi}
  return()
endif()

# Transparent deps of THIS flavor (other flavors in the same prefix may
# have different ones), then the flavor's targets. The deps file runs
# nested find_dependency() calls whose Configs come from this same
# template and (re)use the same _FC_* variables in this scope — stash
# the targets path in a package-unique variable across that include.
set(_FC_targets_file_${ARG_EXPORT} \"\${_FC_targets_file}\")
include(\"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}-deps-${_consumer_tag_expr}.cmake\")
include(\"\${_FC_targets_file_${ARG_EXPORT}}\")
unset(_FC_targets_file_${ARG_EXPORT})

unset(_FC_consumer_family)
unset(_FC_consumer_tag)
unset(_FC_targets_file)
${_cleanup_mpi}
")
  else()
    set(_config_content "\
# ${_config_name}Config.cmake
# Auto-generated by FortranCompiler.cmake
#
# This package's archive contains no Fortran objects and is not bound
# to an MPI ABI, so it is compatible across Fortran compilers and MPI
# flavors — a single untagged targets file suffices and no
# consumer-side tag detection is performed.

cmake_minimum_required(VERSION 3.12)

include(CMakeFindDependencyMacro)
${_mpi_dep_block}${_deps_block}
include(\"\${CMAKE_CURRENT_LIST_DIR}/${ARG_EXPORT}.cmake\")
")
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
