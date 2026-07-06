# eplinalg — library-assembly helper functions
#
# include()'d from the staged top-level CMakeLists.txt (same directory
# scope — variables and targets defined here behave exactly as if the
# content were inline). Staged flat next to CMakeLists.txt by
# ``migrator stage`` (see the cmake-glue copy lists in
# ``src/migrator/staging.py``).

# Relax the Fortran line-length limits on <target>. The migrator can
# lengthen tokens (e.g. MPI_DOUBLE_COMPLEX → MPI_C_LONG_DOUBLE_COMPLEX
# for kind10), pushing fixed-form lines past column 72 and free-form
# lines past column 132; disable both limits so migrated source
# compiles regardless of growth. Also applied to the pristine genuine
# solvers to stay robust to compiler defaults.
function(fortran_relax_line_length target)
    if(FORTRAN_COMPILER_FAMILY STREQUAL "gfortran"
       OR FORTRAN_COMPILER_FAMILY STREQUAL "flang")
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>
            $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>)
    elseif(FORTRAN_COMPILER_FAMILY STREQUAL "intel")
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:Fortran>:-extend-source>)
    endif()
endfunction()

# Compile <target>'s C sources as C++17 (multifloats mode: the migrated
# C uses operator overloading on the double-double types). Any source
# files passed after the target are re-marked LANGUAGE CXX; the target
# gets the permissive flags the K&R-era code needs, plus the macros
# that skip the C++ MPI bindings — mpicxx.h defines template classes
# that conflict with the ``extern "C"`` wrap the c_migrator post-pass
# injects around .c bodies (the C bindings in mpi.h itself are all we
# need).
function(compile_c_as_cxx target)
    if(ARGN)
        set_source_files_properties(${ARGN} PROPERTIES LANGUAGE CXX)
    endif()
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_compile_options(${target} PRIVATE -fpermissive -Wno-write-strings)
    target_compile_definitions(${target} PRIVATE
        MPICH_SKIP_MPICXX OMPI_SKIP_MPICXX)
endfunction()

# target_link_libraries(<target> PUBLIC <dep>...) where both the target
# and each dep are linked only if they exist in this configuration —
# the per-library sections use this to pair a migrated archive with its
# standard-precision counterpart and its staged dependencies without
# repeating the TARGET guards.
function(link_if_present target)
    if(NOT TARGET ${target})
        return()
    endif()
    foreach(_dep ${ARGN})
        if(TARGET ${_dep})
            target_link_libraries(${target} PUBLIC ${_dep})
        endif()
    endforeach()
endfunction()

# Build a stock, standard-precision archive (S/D/C/Z entry points) from
# vendored upstream sources. The target name is the bare library name
# (``blas``, ``lapack``, ...) — no ``_std`` suffix — so consumers can
# write ``target_link_libraries(... blas)`` and get the upstream
# DGEMM_/ZGEMM_/etc. symbols. The migrated archive (``${LIB_PAIR_PREFIX}<name>``)
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

    # The precision target carries the full family-pair prefix (ey/qx/mw)
    # so the target name, the emitted archive filename and the installed
    # package name all agree: target eylapack → libeylapack.a →
    # find_package(eylapack). (LIB_PREFIX / LIB_PREFIX_COMPLEX remain the
    # single-letter SYMBOL prefixes — EGEMM/YGEMM etc. — and are not used
    # for target names.)
    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PAIR_PREFIX}${lib_name}")

    if(_common)
        add_library(${_common_target} STATIC ${_common})
        fortran_module_layout(${_common_target})
    endif()

    add_library(${_precision_target} STATIC ${_precision})
    fortran_module_layout(${_precision_target})
    if(TARGET ${_common_target})
        target_link_libraries(${_precision_target} PUBLIC ${_common_target})
    endif()

    fortran_relax_line_length(${_precision_target})
    if(TARGET ${_common_target})
        fortran_relax_line_length(${_common_target})
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

    # Family-pair-prefixed target name, matching the archive filename and
    # package name (see add_migrated_fortran_library).
    set(_common_target "${lib_name}_common")
    set(_precision_target "${LIB_PAIR_PREFIX}${lib_name}")

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
            compile_c_as_cxx(${_dual_target} ${_dsrc})
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
        compile_c_as_cxx(${_precision_target} ${_all_sources})
        if(TARGET ${_common_target})
            compile_c_as_cxx(${_common_target})
        endif()
        # Link multifloats + bridge for C libraries in multifloats mode.
        # ``multifloats`` is the header-only C++ double-double core: it is
        # only needed to *compile* this target's own C++ sources — the
        # instantiated symbols are baked into the resulting archive — so it
        # is a BUILD_INTERFACE-only dependency. Leaking ``multifloats::``
        # into the install export would make find_package(${_precision_target})
        # fail, because the FetchContent'd core is not itself an installed,
        # find_package-able package.
        if(TARGET multifloats)
            target_link_libraries(${_precision_target} PUBLIC $<BUILD_INTERFACE:multifloats>)
        endif()
        if(TARGET multifloats_mpi)
            target_link_libraries(${_precision_target} PUBLIC multifloats_mpi)
            if(TARGET ${_common_target})
                target_link_libraries(${_common_target} PUBLIC multifloats_mpi)
            endif()
        endif()
    endif()
endfunction()
