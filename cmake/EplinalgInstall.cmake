# eplinalg — install/export machinery and rules
#
# include()'d from the staged top-level CMakeLists.txt (same directory
# scope — variables and targets defined here behave exactly as if the
# content were inline). Staged flat next to CMakeLists.txt by
# ``migrator stage`` (see the cmake-glue copy lists in
# ``src/migrator/staging.py``).

# ── Install ──────────────────────────────────────────────────────────
# Each precision library gets its own export set and Config.cmake,
# so consumers can do: find_package(qxblas), find_package(qxlapack), etc.
# The common library is bundled into the same export as its precision library.

# Libraries whose *objects* are MPI-ABI-bound (reference MPI symbols or
# bake mpi.h constants) — their installs pass ``MPI`` so the archive
# filename picks up a ``-${MPI_TAG}`` suffix and OpenMPI/MPICH/Intel-MPI
# variants can coexist in the same install prefix. Per an nm audit of
# the installed archives: only blacs and mumps objects reference MPI;
# pblas/pbblas/scalapack reach MPI exclusively through BLACS (their
# archives are byte-identical across MPI flavors), so they carry no
# flavor tag. (ptscotch is also MPI-ABI-tied and carries -${MPI_TAG},
# but it is plain-installed alongside scotch in the mumps block below —
# where the tag is baked into its OUTPUT_NAME at target creation — not
# through this list.)
set(_MPI_ABI_LIBS blacs mumps)

# Libraries that cannot build without MPI at all (their sources include
# mpi.h or call BLACS) — installs are skipped entirely when MPI is
# absent, independent of whether the resulting archive is ABI-bound.
set(_MPI_REQUIRED_LIBS blacs pblas pbblas scalapack mumps)

# _install_target_and_modules(<target> <fortran_install_library args…>):
# every install site pairs fortran_install_library with the module
# install (Fortran .mod files ship whenever the target compiled any).
# ARGN (EXPORT/NAMESPACE/MPI/DEPENDS) is forwarded verbatim.
function(_install_target_and_modules target)
    fortran_install_library(${target} ${ARGN})
    get_target_property(_moddir ${target} Fortran_MODULE_DIRECTORY)
    if(_moddir)
        fortran_install_modules(${target})
    endif()
endfunction()

# _eplinalg_pkg_name(<kind> <lib_name> <out_var>): package name a shared
# archive is installed under — _eplinalg_pkg_name(Std blas V) →
# eplinalgStdBlas, _eplinalg_pkg_name(Common mumps V) →
# eplinalgCommonMumps. The capitalization keeps the bare lib name
# visually intact while satisfying CMake's package-name conventions
# (the leading capital aids readability when grep'ing for
# find_package(eplinalg...) calls).
function(_eplinalg_pkg_name kind lib_name out_var)
    string(SUBSTRING "${lib_name}" 0 1 _first)
    string(SUBSTRING "${lib_name}" 1 -1 _rest)
    string(TOUPPER "${_first}" _first)
    set(${out_var} "eplinalg${kind}${_first}${_rest}" PARENT_SCOPE)
endfunction()

# _return_if_already_installed(<key>): once-per-configure latch for the
# shared-package installers — several precision libs can trigger the
# same shared install in one configure; only the first wins. A macro,
# so the return() unwinds the CALLING function.
macro(_return_if_already_installed key)
    get_property(_done_once GLOBAL PROPERTY ${key})
    if(_done_once)
        return()
    endif()
    set_property(GLOBAL PROPERTY ${key} TRUE)
endmacro()

# _append_linked_bridge_deps(<target> <out_list>): append whichever
# first-party MPI bridge packages (multifloats_mpi(_f) for m/w,
# quad_mpi(_f) for q/x — each installed as its own Config package in
# the NEEDS_* blocks at the end of this file) <target> actually links.
# A generated Targets file that references eplinalg::<bridge> without
# the matching find_dependency leaves the imported target undefined at
# consumer time. Reading the target's real LINK_LIBRARIES beats
# re-deriving the per-site NEEDS_* guards.
function(_append_linked_bridge_deps target out_list)
    set(_deps ${${out_list}})
    get_target_property(_ll ${target} LINK_LIBRARIES)
    if(_ll)
        foreach(_b multifloats_mpi multifloats_mpi_f quad_mpi quad_mpi_f)
            if(TARGET ${_b} AND "${_b}" IN_LIST _ll)
                list(APPEND _deps ${_b})
            endif()
        endforeach()
    endif()
    set(${out_list} "${_deps}" PARENT_SCOPE)
endfunction()

# install_library_headers(<lib_name>): install any ``*.h`` headers
# under ``${lib_name}/src/`` to ``include/${lib_name}/`` and extend
# every related target's INTERFACE include directories so consumers
# linking against the installed library pick up the include path via
# ``find_package``. Per-library subdirs are kept as a defensive
# convention against header-name clashes across libraries (PBLAS and
# the ScaLAPACK C layer both carry an internal ``pblas.h`` with
# different contents, but the latter's internal headers are never
# installed — no two installed libraries share a header basename).
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
    foreach(_t ${LIB_PAIR_PREFIX}${lib_name} ${lib_name} ${lib_name}_common)
        if(TARGET ${_t})
            target_include_directories(${_t}
                INTERFACE $<INSTALL_INTERFACE:include/${lib_name}>)
        endif()
    endforeach()
    install(DIRECTORY ${_src}/
        DESTINATION include/${lib_name}
        FILES_MATCHING PATTERN "*.h")
endfunction()

# Install the standard-precision archive ``lib_name`` (e.g. ``blas``,
# ``lapack``) into its own package ``eplinalgStd<LibName>``. The archive
# is byte-identical across per-precision build invocations (vanilla
# Netlib sources, no real-promotion flag, no precision-specific
# defines), so re-installing on top of an existing prefix overwrites
# matching content — which is exactly what consumers want: one
# ``eplinalg::blas`` target shared by ``find_package(qxblas)``,
# ``find_package(eyblas)``, ``find_package(mwblas)``, ... avoiding the
# multi-target redefinition that bundling the std archive in each
# per-precision export used to cause.
function(_install_standard_archive lib_name)
    if(NOT TARGET ${lib_name})
        return()
    endif()
    _eplinalg_pkg_name(Std ${lib_name} _pkg)
    set(_export "${_pkg}Targets")
    _return_if_already_installed(_FC_STD_INSTALLED_${_pkg})

    set(_mpi_arg "")
    if("${lib_name}" IN_LIST _MPI_ABI_LIBS)
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
                _eplinalg_pkg_name(Std "${_d}" _dep_pkg)
                if(NOT "${_dep_pkg}" STREQUAL "${_pkg}")
                    list(APPEND _std_deps "${_dep_pkg}")
                endif()
            endif()
        endforeach()
    endif()

    _install_target_and_modules(${lib_name}
        EXPORT ${_export}
        NAMESPACE eplinalg::
        ${_mpi_arg}
        DEPENDS ${_std_deps})
endfunction()

# ── Shared arith-agnostic packages ───────────────────────────────────
# The full 10-arithmetic release is built from three separately-migrated
# trees (kind10→e/y, kind16→q/x, multifloats→m/w) installed into one
# prefix. Every tree ALSO rebuilds the identical arith-agnostic shared
# targets (the ``*_common`` archives, the genuine s/c/d/z solvers, the
# ordering libraries). If each tree bundles those into its OWN typed
# export set, the generated ``*Targets.cmake`` collide at find_package
# time ("Some (but not all) targets in this export set were already
# defined") and the release cannot be find_package'd for all ten
# precisions at once. So — mirroring ``_install_standard_archive``
# for the Netlib archives — each shared target is exported ONCE, in its
# own package, and the precision Configs ``find_dependency`` it. The
# content-driven archive filenames make the three trees emit
# byte-identical ``.a`` files, so the shared package is a genuine
# last-writer-wins no-op on re-install.

# Install ``${lib_name}_common`` into its own eplinalgCommon<Lib> package
# (once per configure). MPI-tagged for the MPI-dependent libraries.
# mumps_common PUBLIC-links the ordering leaves, so its Config
# find_dependency's the ordering package.
function(_install_shared_common lib_name)
    set(_ct "${lib_name}_common")
    if(NOT TARGET ${_ct})
        return()
    endif()
    _eplinalg_pkg_name(Common ${lib_name} _pkg)
    _return_if_already_installed(_FC_COMMON_INSTALLED_${_pkg})

    set(_mpi_arg "")
    if("${lib_name}" IN_LIST _MPI_ABI_LIBS)
        set(_mpi_arg "MPI")
    endif()
    set(_deps "")
    if("${lib_name}" STREQUAL "mumps" AND
       (TARGET pord OR TARGET scotch OR TARGET metis))
        list(APPEND _deps eplinalgOrdering)
    endif()

    # C ``_common`` archives built in multifloats mode PUBLIC-link the
    # multifloats_mpi bridge (see add_migrated_c_library's dual-interface
    # block); mumps_common likewise reaches the first-party bridges.
    _append_linked_bridge_deps(${_ct} _deps)

    _install_target_and_modules(${_ct}
        EXPORT ${_pkg}Targets
        NAMESPACE eplinalg::
        ${_mpi_arg}
        DEPENDS ${_deps})
endfunction()

# Install the MUMPS ordering leaf archives (pord/metis/scotch/esmumps
# + the MPI-tied ptscotch) into a single
# eplinalgOrdering package (once per configure). These are pure-C leaf
# archives whose OUTPUT_NAMEs are baked at target creation
# (``libpord_mumps.a`` …) — they must NOT go through
# fortran_install_library (which would re-tag the filename), so the
# export + a minimal Config are emitted by hand. mumps_common and the
# genuine/typed solvers reach these via find_dependency(eplinalgOrdering).
function(_install_ordering_package)
    if(NOT TARGET pord AND NOT TARGET scotch AND NOT TARGET metis)
        return()
    endif()
    _return_if_already_installed(_FC_ORDERING_INSTALLED)

    set(_ord "")
    foreach(_t pord metis scotch esmumps ptscotch)
        if(TARGET ${_t})
            list(APPEND _ord ${_t})
        endif()
    endforeach()
    install(TARGETS ${_ord} EXPORT eplinalgOrderingTargets)
    install(EXPORT eplinalgOrderingTargets
        FILE eplinalgOrderingTargets.cmake
        NAMESPACE eplinalg::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/eplinalgOrdering)

    # ptscotch links MPI::MPI_C; only then does the Config need MPI.
    set(_ord_mpi_dep "")
    if(TARGET ptscotch)
        set(_ord_mpi_dep "find_dependency(MPI COMPONENTS C)")
    endif()
    set(_cfg "\
# eplinalgOrderingConfig.cmake — auto-generated by cmake/CMakeLists.txt
# MUMPS ordering leaf archives (PORD / METIS / Scotch / [PT-Scotch]).
# Pure-C archives with a stable ABI, so no compiler-tag matching is
# needed — a single untagged targets file suffices.
cmake_minimum_required(VERSION 3.12)
include(CMakeFindDependencyMacro)
${_ord_mpi_dep}
include(\"\${CMAKE_CURRENT_LIST_DIR}/eplinalgOrderingTargets.cmake\")
")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/eplinalgOrderingConfig.cmake "${_cfg}")
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/eplinalgOrderingConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/eplinalgOrdering)
endfunction()

# Install the genuine single+double MUMPS solvers (dzmumps/scmumps +
# their RESCAN umbrellas) into eplinalgGenuineMumps (once per configure).
# They PUBLIC-link the shared mumps_common → find_dependency(eplinalgCommonMumps).
function(_install_genuine_mumps)
    if(NOT TARGET dzmumps AND NOT TARGET scmumps)
        return()
    endif()
    _return_if_already_installed(_FC_GENUINE_INSTALLED)

    foreach(_g dzmumps scmumps)
        if(TARGET ${_g})
            _install_target_and_modules(${_g}
                EXPORT eplinalgGenuineMumpsTargets
                NAMESPACE eplinalg:: MPI
                DEPENDS eplinalgCommonMumps)
        endif()
    endforeach()
    foreach(_gf dzmumps_full scmumps_full)
        if(TARGET ${_gf})
            install(TARGETS ${_gf} EXPORT eplinalgGenuineMumpsTargets)
        endif()
    endforeach()
endfunction()

function(_install_library_pair lib_name)
    set(_precision_target "${LIB_PAIR_PREFIX}${lib_name}")
    set(_common_target "${lib_name}_common")
    set(_export "${_precision_target}Targets")

    if(NOT TARGET ${_precision_target})
        return()
    endif()

    set(_install_mpi_arg "")
    if("${lib_name}" IN_LIST _MPI_ABI_LIBS)
        set(_install_mpi_arg "MPI")
    endif()
    if("${lib_name}" IN_LIST _MPI_REQUIRED_LIBS)
        # MPI-requiring libs only compile when MPI is available; the
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

    # Factor the shared arith-agnostic targets into their own packages
    # BEFORE the precision export is written, so its INTERFACE
    # references (eplinalg::${lib}_common, eplinalg::pord, …) resolve as
    # cross-export references into those standalone packages rather than
    # being bundled into — and colliding across — the typed export sets.
    if("${lib_name}" STREQUAL "mumps")
        _install_ordering_package()
    endif()
    _install_shared_common(${lib_name})
    if("${lib_name}" STREQUAL "mumps")
        _install_genuine_mumps()
    endif()

    # If the std archive exists, its package is a transparent dep of
    # the precision Config so consumers only need find_package(qxblas)
    # and ``eplinalg::blas`` resolves automatically.
    set(_precision_deps "")
    if(TARGET ${lib_name} AND NOT "${lib_name}" STREQUAL "${_precision_target}")
        _eplinalg_pkg_name(Std ${lib_name} _std_pkg)
        list(APPEND _precision_deps "${_std_pkg}")
    endif()

    # A subset of precision targets (the blacs / scalapack / mumps
    # sections) PUBLIC-link the first-party MPI bridges —
    # multifloats_mpi(_f) for the m/w double-double MPI ops, quad_mpi(_f)
    # for the q/x __float128 ops. Each bridge is installed as its own
    # Config package (see the NEEDS_MULTIFLOATS / NEEDS_QUAD_MPI install
    # blocks below). Whichever bridge the precision target actually links
    # must be a find_dependency of its Config, else the generated typed
    # Targets file carries an ``eplinalg::multifloats_mpi`` (etc.)
    # INTERFACE_LINK_LIBRARIES entry that resolves to an undefined
    # imported target when a consumer links the package. Read the target's
    # real LINK_LIBRARIES (bare target names) rather than re-deriving the
    # per-site NEEDS_* guards.
    _append_linked_bridge_deps(${_precision_target} _precision_deps)

    # The precision archive PUBLIC-links its ``_common`` (factored into
    # eplinalgCommon<Lib>), so its Config find_dependency's it.
    if(TARGET ${_common_target})
        _eplinalg_pkg_name(Common ${lib_name} _common_pkg)
        list(APPEND _precision_deps ${_common_pkg})
    endif()

    # MUMPS precision archives additionally reach the ordering leaves
    # (${LIB_PAIR_PREFIX}mumps PUBLIC-links ptscotch directly, and
    # mumps_common the rest), and — for consumer convenience — the
    # genuine s/c/d/z solver package.
    if("${lib_name}" STREQUAL "mumps")
        if(TARGET pord OR TARGET scotch OR TARGET metis)
            list(APPEND _precision_deps eplinalgOrdering)
        endif()
        if(TARGET dzmumps OR TARGET scmumps)
            list(APPEND _precision_deps eplinalgGenuineMumps)
        endif()
    endif()

    # Install the precision-specific archive. fortran_install_library
    # writes ${_export}'s Config.cmake on this first call, emitting a
    # find_dependency() for each transparent package in _precision_deps.
    _install_target_and_modules(${_precision_target}
        EXPORT ${_export}
        NAMESPACE eplinalg::
        ${_install_mpi_arg}
        DEPENDS ${_precision_deps})

    # The MUMPS typed umbrella (${LIB_PAIR_PREFIX}mumps_full — a RESCAN
    # wrapper over the typed solver + shared mumps_common) is
    # precision-specific, so it stays in the typed export. The shared C
    # runtime (mumps_common), the ordering leaves and the genuine solvers
    # were already factored into their standalone packages above.
    if("${lib_name}" STREQUAL "mumps" AND TARGET ${LIB_PAIR_PREFIX}mumps_full)
        install(TARGETS ${LIB_PAIR_PREFIX}mumps_full EXPORT ${_export})
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
if(NOT BASELINE_BUILD AND TARGET ${LIB_PAIR_PREFIX}mumps)
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
    # The double-double *Fortran* module ``multifloatsf`` (module name
    # ``multifloats``, symbols ``__multifloats_MOD_*``) is a genuine link
    # dependency of every m/w Fortran archive: those archives leave the
    # module procedures as *undefined* references (we link it BUILD_INTERFACE
    # so they compile against its .mod, but the symbols are NOT baked in the
    # way the header-only C++ core is). A find_package consumer therefore has
    # to link libmultifloatsf itself, so it must be installed as its own
    # package. The FetchContent'd multifloats sub-build gates its own install
    # off (MULTIFLOATSF_INSTALL_PRECOMPILED_MOD, and it is EXCLUDE_FROM_ALL),
    # so we install the target from here. Non-MPI, clean interface (PRIVATE
    # includes only), so no DEPENDS / MPI tag. Built on demand as a link
    # dependency of the m/w targets (which ARE in ALL), so the archive exists
    # by install time despite EXCLUDE_FROM_ALL.
    # The C++ double-double *core* ``multifloats`` (extern "C" adddd / cadddd /
    # sindd / csqrtdd … from multifloats_math.cc) is a genuine link dependency:
    # the Fortran module procedures in libmultifloatsf call straight into these
    # C-ABI functions, and they are NOT baked into any m/w archive (only
    # header-instantiated template ops are). The FetchContent'd sub-build already
    # installs it as package ``multifloats`` under namespace ``multifloats::``
    # (add_subdirectory EXCLUDE_FROM_ALL excludes it from the default *build*
    # target, but the sub-build's install() rules still run), so we do NOT
    # re-install it here — doing so would place ``multifloats`` in a second
    # export set and break install(EXPORT multifloatsfTargets), which requires
    # its dependency to live in exactly one. We only wire multifloatsf's Config
    # to find it (below).

    if(TARGET multifloatsf)
        # multifloatsf's exported interface references multifloats::multifloats
        # (PUBLIC-linked C++ core), so emit find_dependency(multifloats) into its
        # Config — otherwise find_package(multifloatsf) fails with "imported
        # targets referenced, but missing: multifloats::multifloats".
        _install_target_and_modules(multifloatsf EXPORT multifloatsfTargets NAMESPACE eplinalg:: DEPENDS multifloats)
    endif()

    # multifloats_mpi's C objects register the custom MPI ops and are
    # MPI-ABI-bound → flavor tag. multifloats_mpi_f is a single Fortran
    # binding shim with no MPI references of its own → compiler tag
    # only (auto-detected); it PUBLIC-links multifloats_mpi, so its
    # Targets file references eplinalg::multifloats_mpi — emit the
    # matching find_dependency into its Config so a consumer that loads
    # only multifloats_mpi_f still resolves the base bridge.
    if(TARGET multifloats_mpi)
        _install_target_and_modules(multifloats_mpi EXPORT multifloats_mpiTargets NAMESPACE eplinalg:: MPI)
    endif()
    if(TARGET multifloats_mpi_f)
        _install_target_and_modules(multifloats_mpi_f EXPORT multifloats_mpi_fTargets NAMESPACE eplinalg:: DEPENDS multifloats_mpi)
    endif()
endif()

if(NEEDS_QUAD_MPI)
    # quad_mpi's C objects register the custom MPI ops and are
    # MPI-ABI-bound → flavor tag. quad_mpi_f is a single Fortran binding
    # shim with no MPI references of its own → compiler tag only
    # (auto-detected); it PUBLIC-links quad_mpi, so its Targets file
    # references eplinalg::quad_mpi — emit find_dependency(quad_mpi)
    # into its Config so loading only quad_mpi_f still resolves the
    # base bridge.
    if(TARGET quad_mpi)
        _install_target_and_modules(quad_mpi EXPORT quad_mpiTargets NAMESPACE eplinalg:: MPI)
    endif()
    if(TARGET quad_mpi_f)
        _install_target_and_modules(quad_mpi_f EXPORT quad_mpi_fTargets NAMESPACE eplinalg:: DEPENDS quad_mpi)
    endif()
endif()
