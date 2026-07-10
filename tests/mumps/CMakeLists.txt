# Differential precision tests for migrated MUMPS.
#
# Mirrors tests/lapack/CMakeLists.txt with two extras:
#   1. C-side bridge — upstream mumps_c.c is compiled twice (Q + X arith)
#      against target_${TARGET_NAME}/c/include/'s mumps_c_types_extended.h
#      (extends upstream mumps_c_types.h to the target widths), plus
#      type-agnostic upstream helpers (mumps_common.c, mumps_addr.c, IO
#      / save / thread / metis / pord / scotch). See tests/mumps/TODO.md
#      B2 for the override design.
#   2. Every test runs via ``mpiexec -n 1`` — the migrated qxmumps Fortran
#      calls MPI primitives (mpi_init_, mpi_bcast_, ...) regardless of
#      NPROCS. Real MPI is required at runtime; libmpiseq exists as an
#      alternative (B3 in TODO.md) but its full sequential-link path is
#      not yet wired in.
#
# Outputs: per-test JSON reports under
# ``${CMAKE_BINARY_DIR}/precision_reports/<test_name>.<target>.json``.

# ── Gate checks ──────────────────────────────────────────────────────
if(NOT TARGET ${LIB_PAIR_PREFIX}mumps)
    message(STATUS "tests/mumps: skipping — ${LIB_PAIR_PREFIX}mumps not built")
    return()
endif()
if(NOT HAVE_REAL16)
    message(STATUS "tests/mumps: skipping — HAVE_REAL16 required")
    return()
endif()
if(NOT TARGET reflapack_quad)
    message(STATUS "tests/mumps: skipping — reflapack_quad not available "
                   "(tests/lapack must be enabled)")
    return()
endif()
if(NOT MPI_Fortran_FOUND OR NOT MPI_C_FOUND)
    message(STATUS "tests/mumps: skipping — MPI Fortran+C not found")
    return()
endif()

set(_wrapper_dir ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME})
if(NOT IS_DIRECTORY ${_wrapper_dir})
    message(STATUS "tests/mumps: no wrapper for target '${TARGET_NAME}' "
                   "at ${_wrapper_dir} — skipping")
    return()
endif()

set(_MUMPS_C_SRC ${PROJECT_SOURCE_DIR}/_mumps_upstream_src)
set(_MUMPS_C_INC ${PROJECT_SOURCE_DIR}/_mumps_upstream_include)
if(NOT IS_DIRECTORY ${_MUMPS_C_SRC} OR NOT IS_DIRECTORY ${_MUMPS_C_INC})
    message(STATUS "tests/mumps: skipping — upstream MUMPS C sources not "
                   "staged at ${_MUMPS_C_SRC} / ${_MUMPS_C_INC}")
    return()
endif()

# Per-target bridge include directory carries two independent artifacts:
# mumps_c_types_extended.h — force-included ONLY into the mumps_c.c bridge
# objects to override upstream's ``double`` widths — and the standalone
# consumer headers ${LIB_PREFIX}mumps_c.h / ${LIB_PREFIX_COMPLEX}mumps_c.h,
# which carry their own widened widths inline and do NOT include the
# extended header. ${PFX} / ${PFX_C} are the upper-case arithmetic
# prefixes used in the bridge's macro renames. The target-agnostic
# mumps_int_def.h sits in the shared c/include/.
set(_BRIDGE_INC_TARGET ${CMAKE_CURRENT_SOURCE_DIR}/target_${TARGET_NAME}/c/include)
set(_BRIDGE_INC_SHARED ${CMAKE_CURRENT_SOURCE_DIR}/c/include)
if(NOT IS_DIRECTORY ${_BRIDGE_INC_TARGET})
    message(STATUS "tests/mumps: skipping — no target_${TARGET_NAME}/c/include")
    return()
endif()
string(TOUPPER "${LIB_PREFIX}" PFX)
string(TOUPPER "${LIB_PREFIX_COMPLEX}" PFX_C)

# ── Output directory for JSON reports ────────────────────────────────
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)

# ── Dedicated module directory ───────────────────────────────────────
# Same prec_kinds / compare / prec_report / test_data module names as
# tests/blas and tests/lapack, so isolate the .mod outputs to avoid
# collision when multiple test trees build together.
set(_mp_moddir "${PROJECT_BINARY_DIR}/fmod_tests/${FORTRAN_MOD_COMPAT_TAG}/mumps_tests")
file(MAKE_DIRECTORY ${_mp_moddir})

# ── C-side bridge ─────────────────────────────────────────────────────
# Per-arithmetic dispatch — compile mumps_c.c twice (real + complex)
# as OBJECT libraries so they can be folded into the single combined
# archive below. Override Fortran symbols so ${LIB_PREFIX}mumps_c()
# lands in the bridge but the dmumps_assign_*/dmumps_set_tmp_ptr_c_
# helpers stay under their D-prefix (because the migrated
# ${LIB_PREFIX}mumps_f77_ Fortran expects them that way — same
# convention for Z/${LIB_PREFIX_COMPLEX}).
add_library(${LIB_PREFIX}mumps_c_bridge OBJECT ${_MUMPS_C_SRC}/mumps_c.c)
target_include_directories(${LIB_PREFIX}mumps_c_bridge PRIVATE
    ${_BRIDGE_INC_TARGET} ${_BRIDGE_INC_SHARED} ${_MUMPS_C_INC} ${_MUMPS_C_SRC}
)
target_compile_definitions(${LIB_PREFIX}mumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_d Add_
    dmumps_f77_=${LIB_PREFIX}mumps_f77_
    dmumps_set_tmp_ptr_=${LIB_PREFIX}mumps_set_tmp_ptr_
    dmumps_c=${LIB_PREFIX}mumps_c
    DMUMPS_STRUC_C=${PFX}MUMPS_STRUC_C
)
# Force-include the per-target mumps_c_types_extended.h BEFORE
# mumps_c.c pulls in mumps_c_types.h via mumps_common.h — the extended
# header includes upstream's (guarded) then overrides its ``double``
# widths, so the whole translation unit sees the target types.
target_compile_options(${LIB_PREFIX}mumps_c_bridge PRIVATE
    -include mumps_c_types_extended.h)

add_library(${LIB_PREFIX_COMPLEX}mumps_c_bridge OBJECT ${_MUMPS_C_SRC}/mumps_c.c)
target_include_directories(${LIB_PREFIX_COMPLEX}mumps_c_bridge PRIVATE
    ${_BRIDGE_INC_TARGET} ${_BRIDGE_INC_SHARED} ${_MUMPS_C_INC} ${_MUMPS_C_SRC}
)
target_compile_definitions(${LIB_PREFIX_COMPLEX}mumps_c_bridge PRIVATE
    MUMPS_ARITH=MUMPS_ARITH_z Add_
    zmumps_f77_=${LIB_PREFIX_COMPLEX}mumps_f77_
    zmumps_set_tmp_ptr_=${LIB_PREFIX_COMPLEX}mumps_set_tmp_ptr_
    zmumps_c=${LIB_PREFIX_COMPLEX}mumps_c
    ZMUMPS_STRUC_C=${PFX_C}MUMPS_STRUC_C
)
target_compile_options(${LIB_PREFIX_COMPLEX}mumps_c_bridge PRIVATE
    -include mumps_c_types_extended.h)

# mumps_addr.c stays at upstream double precision because the migrated
# Fortran callers (qtools.F's QMUMPS_UPDATE_PROGRESS, xtools.F's
# XMUMPS_UPDATE_PROGRESS) pass DOUBLE PRECISION OPELI — the keep-kind
# manifest preserved DP for those callers. Compiling mumps_addr.c with
# the extended types would make MUMPS_SETDVAL_ADDR_C write the
# wider value into the 8-byte caller buffer (silent corruption →
# SIGSEGV inside ${LIB_PREFIX}mumps_update_progress_).
add_library(mumps_c_addr_obj OBJECT ${_MUMPS_C_SRC}/mumps_addr.c)
# Include path order matters: list _MUMPS_C_INC FIRST so upstream's
# `double`-typed mumps_c_types.h gets resolved by mumps_addr.c's
# `#include "mumps_c_types.h"` chain (the target dir with the extended
# header is not on this object's path at all). _BRIDGE_INC_SHARED is still on
# the search path so upstream's `#include "mumps_int_def.h"` falls
# through to our generated stand-in (upstream ships only a .in
# template).
target_include_directories(mumps_c_addr_obj PRIVATE
    ${_MUMPS_C_INC} ${_MUMPS_C_SRC} ${_BRIDGE_INC_SHARED})
target_compile_definitions(mumps_c_addr_obj PRIVATE Add_)

# Combined C-side static archive. Bundles the per-arithmetic bridge
# objects with the type-agnostic upstream runtime (mumps_common.c,
# IO, save/restore, thread, numa, flytes, config, metis, pord, scotch)
# into ONE archive. Doing this rather than two archives sidesteps a
# static-link ordering problem: the migrated qxmumps Fortran (built
# upstream of tests/) calls into both the bridge's
# dmumps_set_tmp_ptr_c_ and the runtime's mumps_size_c_ /
# mumps_low_level_*, and a separate runtime archive listed before
# libqmumps would have its objects skipped at scan time and not
# re-pulled when libqmumps's later objects request them.
add_library(mumps_c_bridge STATIC
    $<TARGET_OBJECTS:${LIB_PREFIX}mumps_c_bridge>
    $<TARGET_OBJECTS:${LIB_PREFIX_COMPLEX}mumps_c_bridge>
    $<TARGET_OBJECTS:mumps_c_addr_obj>
    ${CMAKE_CURRENT_SOURCE_DIR}/c/c_parity_helpers.c
    ${_MUMPS_C_SRC}/mumps_common.c
    ${_MUMPS_C_SRC}/mumps_io.c
    ${_MUMPS_C_SRC}/mumps_io_basic.c
    ${_MUMPS_C_SRC}/mumps_io_err.c
    ${_MUMPS_C_SRC}/mumps_io_thread.c
    ${_MUMPS_C_SRC}/mumps_save_restore_C.c
    ${_MUMPS_C_SRC}/mumps_register_thread.c
    ${_MUMPS_C_SRC}/mumps_thread.c
    ${_MUMPS_C_SRC}/mumps_thread_affinity.c
    ${_MUMPS_C_SRC}/mumps_numa.c
    ${_MUMPS_C_SRC}/mumps_flytes.c
    ${_MUMPS_C_SRC}/mumps_config_file_C.c
    ${_MUMPS_C_SRC}/mumps_metis.c
    ${_MUMPS_C_SRC}/mumps_metis64.c
    ${_MUMPS_C_SRC}/mumps_metis_int.c
    ${_MUMPS_C_SRC}/mumps_pord.c
    ${_MUMPS_C_SRC}/mumps_scotch.c
    ${_MUMPS_C_SRC}/mumps_scotch64.c
    ${_MUMPS_C_SRC}/mumps_scotch_int.c
)
target_include_directories(mumps_c_bridge PUBLIC
    $<BUILD_INTERFACE:${_BRIDGE_INC_TARGET}>  # FIRST — per-target extended types
    $<BUILD_INTERFACE:${_BRIDGE_INC_SHARED}>  # then mumps_int_def.h fallback
    $<BUILD_INTERFACE:${_MUMPS_C_INC}>        # then upstream include/
    $<BUILD_INTERFACE:${_MUMPS_C_SRC}>        # then upstream src/
)
target_compile_definitions(mumps_c_bridge PRIVATE Add_
    # c_parity_helpers.c references the per-target header path and
    # entry-point names through these macros so the helper compiles
    # for any target.
    "TARGET_REAL_HEADER=\"${LIB_PREFIX}mumps_c.h\""
    "TARGET_COMPLEX_HEADER=\"${LIB_PREFIX_COMPLEX}mumps_c.h\""
    TARGET_REAL_MUMPS_C=${LIB_PREFIX}mumps_c
    TARGET_COMPLEX_MUMPS_C=${LIB_PREFIX_COMPLEX}mumps_c
    TARGET_REAL_STRUC_C=${PFX}MUMPS_STRUC_C
    TARGET_COMPLEX_STRUC_C=${PFX_C}MUMPS_STRUC_C
)
# Override applies to everything in mumps_c_bridge EXCEPT the
# already-built mumps_c_addr_obj (which we added via $<TARGET_OBJECTS>
# so its compile flags are frozen). The non-_addr files are either
# integer-only (mumps_common.c, IO/save/thread/etc. helpers — the
# override is harmless) or they're the per-arithmetic bridge sources
# that NEED the extended types.
target_compile_options(mumps_c_bridge PRIVATE
    -include mumps_c_types_extended.h)
target_link_libraries(mumps_c_bridge PRIVATE MPI::MPI_C)

# PORD ordering (ICNTL(7)=4). The bridge recompiles mumps_pord.c, whose
# MUMPS_PORDF / MUMPS_PORDF_WND bodies are guarded by ``#if defined(pord)``.
# Without -Dpord here the bridge object is empty, and the migrated
# ana_orderings_wrappers_m.F (in mumps_common) references mumps_pordf_ with
# no definition on the test link line. Mirror the shipped mumps_common
# (folded C runtime) wiring: define pord, expose the PORD headers, and link libpord so
# SPACE_ordering et al. resolve in the final test executables. Gated on the
# pord target so the tree still builds when PORD wasn't staged.
if(TARGET pord)
    target_compile_definitions(mumps_c_bridge PRIVATE pord)
    target_include_directories(mumps_c_bridge PRIVATE
        $<TARGET_PROPERTY:pord,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(mumps_c_bridge PUBLIC pord)
endif()

# METIS ordering (ICNTL(7)=5). Same story as PORD: the bridge recompiles
# mumps_metis.c / mumps_metis64.c / mumps_metis_int.c, whose bodies are
# guarded by ``#if defined(metis)`` and forward to the privately-namespaced
# METIS_MUMPS_* symbols. Without -Dmetis here those objects are empty and the
# migrated ana_orderings_wrappers_m.F references them with no definition on
# the test link line. Mirror the shipped mumps_common (folded C runtime) wiring: define metis,
# expose the METIS public header, and link libmetis so METIS_MUMPS_NodeND et
# al. resolve. Gated on the metis target so the tree still builds when METIS
# wasn't staged.
if(TARGET metis)
    target_compile_definitions(mumps_c_bridge PRIVATE metis)
    target_include_directories(mumps_c_bridge PRIVATE
        $<TARGET_PROPERTY:metis,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(mumps_c_bridge PUBLIC metis)
endif()

# Scotch ordering (ICNTL(7)=3). Same story again: the bridge recompiles
# mumps_scotch.c / mumps_scotch64.c / mumps_scotch_int.c, whose bodies are
# guarded by ``#if defined(scotch)`` and call bare SCOTCH_* / esmumps names.
# Without -Dscotch here those objects are empty and ana_orderings_wrappers_m.F
# references them with no definition. Mirror the shipped mumps_common
# (folded C runtime) wiring: define scotch, expose the Scotch include dir, force-include
# scotch_rename_mumps.h so the bare caller names bind to the _mumps-suffixed
# archive symbols, and link scotch/esmumps. Gated on the scotch
# target so the tree still builds when Scotch wasn't staged.
set(_MUMPS_SCOTCH_INC ${PROJECT_SOURCE_DIR}/_mumps_scotch_include)
if(TARGET scotch)
    target_compile_definitions(mumps_c_bridge PRIVATE scotch)
    target_include_directories(mumps_c_bridge PRIVATE
        $<TARGET_PROPERTY:scotch,INTERFACE_INCLUDE_DIRECTORIES>)
    # SHELL: keeps -include + path together so CMake doesn't dedup it
    # against the -include mumps_c_types_extended.h added above (which
    # would collapse both -include tokens and orphan the second header).
    target_compile_options(mumps_c_bridge PRIVATE
        "SHELL:-include ${_MUMPS_SCOTCH_INC}/scotch_rename_mumps.h")
    target_link_libraries(mumps_c_bridge PUBLIC esmumps scotch)
endif()

# PT-Scotch parallel analysis (ICNTL(28)=2, ICNTL(29)=1). The bridge
# recompiles mumps_scotch.c, whose MUMPS_DGRAPHINIT shim calls bare
# SCOTCH_dgraphInit under ``#if defined(ptscotch)`` and includes
# ptscotch.h. Mirror the shipped mumps_common wiring: define ptscotch,
# force-include scotch_rename_pt_mumps.h so that bare distributed call
# binds to the _mumps-suffixed archive symbol, and link the distributed
# increment. ptscotch PUBLIC-depends on scotch, so CMake emits it ahead
# on the static link line (single-pass ld). Gated on the
# ptscotch target so the tree still builds when PT-Scotch wasn't staged
# (libseq / no-MPI columns).
if(TARGET ptscotch)
    target_compile_definitions(mumps_c_bridge PRIVATE ptscotch)
    target_include_directories(mumps_c_bridge PRIVATE
        $<TARGET_PROPERTY:ptscotch,INTERFACE_INCLUDE_DIRECTORIES>)
    target_compile_options(mumps_c_bridge PRIVATE
        "SHELL:-include ${_MUMPS_SCOTCH_INC}/scotch_rename_pt_mumps.h")
    target_link_libraries(mumps_c_bridge PUBLIC ptscotch)
endif()

# ── Common helper modules ────────────────────────────────────────────
add_library(mumps_test_common STATIC
    ../common/prec_kinds.f90
    ../common/compare.f90
    common/prec_report.f90
    common/test_data_mumps.f90
    common/ref_quad_lapack_solve.f90
)
set_target_properties(mumps_test_common PROPERTIES
    Fortran_MODULE_DIRECTORY "${_mp_moddir}"
)
target_include_directories(mumps_test_common PUBLIC
    $<BUILD_INTERFACE:${_mp_moddir}>
)
# prec_report.f90 does ``use mpi`` (to confine reporting to the host rank
# under mpiexec -n N). It needs only mpi.mod at COMPILE time — the actual
# MPI_COMM_RANK symbol is resolved when the test executable links, against
# MPI::MPI_Fortran (plain variant) or libmpiseq (_seq variant). So pull the
# MPI include path WITHOUT taking an MPI link dependency, which would drag
# real MPI into the _seq executables and clash with mpiseq's stubs.
target_include_directories(mumps_test_common PRIVATE
    $<TARGET_PROPERTY:MPI::MPI_Fortran,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(mumps_test_common PUBLIC reflapack_quad)

# ── Per-target wrapper module (kind16-only today) ────────────────────
add_fypp_source(_mumps_wrapper_src
    ${_wrapper_dir}/target_mumps.fypp
    ${CMAKE_CURRENT_SOURCE_DIR}/common
)
add_library(mumps_test_target STATIC ${_mumps_wrapper_src})
set_target_properties(mumps_test_target PROPERTIES
    Fortran_MODULE_DIRECTORY "${_mp_moddir}"
)
target_include_directories(mumps_test_target PUBLIC
    $<BUILD_INTERFACE:${_mp_moddir}>
)
target_link_libraries(mumps_test_target PUBLIC mumps_test_common target_conv)
if(NEEDS_MULTIFLOATS AND TARGET multifloats_mpi_f)
    # On multifloats, target_qmumps / target_xmumps lazily call
    # ``multifloats_mpi_init`` to register the custom MPI datatypes /
    # ops that MUMPS uses internally (MPI_FLOAT64X2, MPI_DD_SUM, ...).
    target_link_libraries(mumps_test_target PUBLIC
        $<BUILD_INTERFACE:multifloats_mpi_f>)
endif()
if(NEEDS_QUAD_MPI AND TARGET quad_mpi_f)
    # On kind16, target_qmumps / target_xmumps lazily call
    # ``quad_mpi_init`` to register the custom quad reduce ops
    # (MPI_QQ_SUM / MPI_XX_* on MPI_REAL16 / MPI_COMPLEX32) that MUMPS
    # uses internally — Intel MPI's builtin MPI_SUM/MAX/MIN have no
    # 16-byte-real kernel. The wrapper needs quad_mpi_f's .mod file to
    # compile and its handles to link.
    target_link_libraries(mumps_test_target PUBLIC
        $<BUILD_INTERFACE:quad_mpi_f>)
endif()
# qxmumps + scalapack + lapack + blas + mumps_c_bridge are NOT attached
# as PUBLIC dependencies because cmake's PUBLIC propagation would place
# them OUTSIDE the per-test LINK_GROUP wrapper. They're listed
# explicitly inside add_mumps_*_test()'s LINK_GROUP so every archive
# carrying cross-references (qxmumps Fortran ↔ bridge's
# dmumps_set_tmp_ptr_c_ / mumps_size_c_ / mumps_low_level_*; bridge ↔
# qmumps_f77_) gets re-scanned until everything resolves. The wrapper
# Fortran source still needs qxmumps's .mod files at compile time, so
# pull just the include dir without taking a link dependency.
get_target_property(_qmumps_moddir ${LIB_PAIR_PREFIX}mumps Fortran_MODULE_DIRECTORY)
if(_qmumps_moddir)
    target_include_directories(mumps_test_target PUBLIC
        $<BUILD_INTERFACE:${_qmumps_moddir}>)
endif()
# Build-order dependency: qmumps_struc_def.mod / xmumps_struc_def.mod
# need to exist before target_mumps.f90 compiles. We didn't set
# qxmumps as a link dependency (PUBLIC propagation would break the
# LINK_GROUP wrapping), so make it explicit here.
add_dependencies(mumps_test_target ${LIB_PAIR_PREFIX}mumps)

# ── Test executables ──────────────────────────────────────────────────
# Rank count for the real-MPI test variant. Default 2 so ctest actually
# exercises the distributed solve path (centralized-solution handling on
# the host, distributed matrix assembly in the dist_input tests) instead
# of only ever running single-rank. The drivers are host-guarded (directly
# or via prec_report), so any np>=1 is correct. The _seq (libmpiseq)
# variant is always single-rank and ignores this. NOTE: run the MPI suite
# with OpenMPI oversubscribe/yield env on core-constrained hosts
# (OMPI_MCA_rmaps_base_oversubscribe=1, OMPI_MCA_mpi_yield_when_idle=1).
set(MUMPS_TEST_NPROCS 2 CACHE STRING
    "MPI rank count for the real-MPI MUMPS test variant")

function(_mumps_register_test name)
    add_test(NAME mumps_${name}
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MUMPS_TEST_NPROCS}
                     ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}>
                     ${MPIEXEC_POSTFLAGS}
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
    set_tests_properties(mumps_${name} PROPERTIES RUN_SERIAL TRUE)
endfunction()

# All static archives that carry cross-references between the migrated
# qxmumps Fortran and the C bridge — collected here so each test wraps
# them in a single LINK_GROUP RESCAN.
set(_mumps_link_group_archives
    mumps_test_target
    mumps_c_bridge
    ${LIB_PAIR_PREFIX}mumps
    # mumps_common carries ana_orderings_wrappers_m, whose only external
    # ordering reference (mumps_pordf_ under -Dpord) is defined in the C
    # bridge. It's a transitive dep of qxmumps and would otherwise land
    # AFTER --end-group, where that reference can't reach back into the
    # bridge. Naming it in the group keeps the mumps_common <-> bridge
    # cycle inside the RESCAN.
    mumps_common
    ${LIB_PAIR_PREFIX}scalapack
    ${LIB_PAIR_PREFIX}pblas
    ${LIB_PAIR_PREFIX}pbblas
    ${LIB_PAIR_PREFIX}ptzblas
    ${LIB_PAIR_PREFIX}blacs
    ${LIB_PAIR_PREFIX}lapack
    ${LIB_PAIR_PREFIX}blas
)
# Under PARALLEL_BLAS, ${LIB_PAIR_PREFIX}blas is an INTERFACE composite.
# $<TARGET_FILE:...> can't resolve interfaces, AND substituting
# _serial into this list creates a CMake target-dep cycle (eylapack
# PUBLIC-links the composite, which links _serial, which is also in
# the group). The fix lives in the per-target foreach below: we keep
# ${LIB_PAIR_PREFIX}blas in the list for add_dependencies(), but in the
# link-line construction we route the blas slot through
# target_link_libraries() (composite) instead of $<TARGET_FILE:>.
# qxxblas is a PUBLIC dep of qxlapack (cmake/CMakeLists.txt:600). qxblas is
# already in the cycle, so qxxblas → qxblas → group → qxxblas closes a
# cycle CMake refuses to build unless qxxblas is named in the LINK_GROUP.
if(TARGET ${LIB_PAIR_PREFIX}xblas)
    list(APPEND _mumps_link_group_archives ${LIB_PAIR_PREFIX}xblas)
endif()
# Ordering archives (ICNTL(7)=4 pord, =5 metis, =3 scotch). The migrated
# solver + bridge + common members reference their symbols (eini/yini_driver.F's
# metis_mumps_*, ana_orderings_wrappers_m's scotchf*_mumps_/metis_mumps_nodend_,
# mumps_metis.c's METIS_MUMPS_*, mumps_pord.c's SPACE_ordering/firstPostorder,
# mumps_scotch.c's SCOTCH_*_mumps), yet mumps_c_bridge/mumps_common only
# PUBLIC-link them — and PUBLIC propagation lands them OUTSIDE the per-test
# LINK_GROUP, where the LINK_GROUP genex drops them. They are NOT named as
# LINK_GROUP targets here: eymumps PUBLIC-links mumps_common which PUBLIC-links
# pord/metis/scotch, so pulling them into the group closes a group<->common
# target-graph cycle. Instead they are appended AFTER the group as raw archive
# files (see _mumps_attach_link_group); collect the target list here, gated on
# their targets, ordered esmumps->scotch for the intra-Scotch chain.
set(_mumps_ordering_archives "")
foreach(_p pord metis esmumps scotch)
    if(TARGET ${_p})
        list(APPEND _mumps_ordering_archives ${_p})
    endif()
endforeach()
# PT-Scotch (distributed, MPI-only) — the -Dptscotch mumps runtime references
# the distributed increment's SCOTCH_dgraph*_mumps / SCOTCHFDGRAPH*_mumps
# symbols. Appended to the same trailing --start-group as the sequential
# ordering archives, so the ptscotch<->scotch cross-references resolve by
# rescan regardless of position. Gated on the ptscotch target: absent in the
# libseq / no-MPI columns, so those links stay sequential-Scotch-only.
if(TARGET ptscotch)
    list(APPEND _mumps_ordering_archives ptscotch)
endif()

# libmpiseq is a complete drop-in for MPI + BLACS + PBLAS + ScaLAPACK
# (single-rank only). The ``_seq`` link variant therefore replaces ALL
# of those real archives with mpiseq alone — keeping mumps + lapack +
# blas + the test harness. No BLACS/scalapack/pblas archives in the
# link, so no symbol collisions with libmpiseq's bundled stubs.
set(_mumps_link_group_archives_seq
    mumps_test_target
    mumps_c_bridge
    ${LIB_PAIR_PREFIX}mumps
    mumps_common   # see _mumps_link_group_archives note (mumps_pordf_)
    ${LIB_PAIR_PREFIX}lapack
    ${LIB_PAIR_PREFIX}blas
)
if(TARGET ${LIB_PAIR_PREFIX}xblas)
    list(APPEND _mumps_link_group_archives_seq ${LIB_PAIR_PREFIX}xblas)
endif()
# (PARALLEL_BLAS workaround: see comment on _mumps_link_group_archives;
# same per-target swap-in-foreach trick applies to _seq.)

# Link-variant registry. Every Fortran/C test source is built twice:
#
#   * ``""``   suffix — links against real MPI (``MPI::MPI_Fortran`` /
#               ``MPI::MPI_C``); test runs under ``mpiexec -n 1``.
#   * ``_seq`` suffix — links against the in-tree ``mpiseq`` archive
#               instead; test runs as a plain binary (no mpiexec).
#
# Both variants pull from the SAME compiled archives — only the MPI
# layer differs. Single-rank numerics should be bit-equivalent across
# the two; the ``_seq`` set exercises the libmpiseq path that
# downstream consumers use when they want a no-MPI executable.
#
# The plain (real-MPI) variant is keyed by the sentinel ``plain`` rather
# than an empty string: a bare "" is a *zero-element* list in CMake, so
# ``foreach(v IN LISTS _mumps_test_variants)`` over "" would iterate zero
# times and silently register NO mumps tests on non-Intel-MPI hosts
# (where ``mpiseq`` is absent). The consumer loops translate ``plain`` ->
# "" for the target-name suffix.
set(_mumps_test_variants "plain")
if(TARGET mpiseq)
    list(APPEND _mumps_test_variants "_seq")
endif()

function(_mumps_apply_link_variant tgt variant want_c)
    if(variant STREQUAL "_seq")
        # mpiseq.a covers both Fortran and C MPI surfaces. Already
        # bundles libmpiseq's own mpi.f / mpic.c plus our supplementary
        # mpiseq_qx_stubs.f / mpiseq_mw_stubs.f90 / mpiseq_c_stubs.c.
        target_link_libraries(${tgt} PRIVATE mpiseq)
        # Test sources still ``use mpi`` / ``#include "mpi.h"`` and
        # reference ``MPI_COMM_WORLD`` / ``MPI_INIT`` etc. that need
        # the real MPI module + headers at COMPILE time. Keep the
        # include path from MPI::MPI_Fortran (and MPI::MPI_C if needed)
        # but DO NOT pull MPI::* into linking — libmpiseq supplies the
        # symbols.
        target_include_directories(${tgt} PRIVATE
            $<TARGET_PROPERTY:MPI::MPI_Fortran,INTERFACE_INCLUDE_DIRECTORIES>)
        if(want_c)
            target_include_directories(${tgt} PRIVATE
                $<TARGET_PROPERTY:MPI::MPI_C,INTERFACE_INCLUDE_DIRECTORIES>)
        endif()
    else()
        target_link_libraries(${tgt} PRIVATE MPI::MPI_Fortran)
        if(want_c)
            target_link_libraries(${tgt} PRIVATE MPI::MPI_C)
        endif()
    endif()
endfunction()

function(_mumps_register_test_variant tgt variant)
    if(variant STREQUAL "_seq")
        # libmpiseq path — plain binary, no mpiexec wrapper.
        add_test(NAME mumps_${tgt}
                 COMMAND $<TARGET_FILE:${tgt}>
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/precision_reports)
        set_tests_properties(mumps_${tgt} PROPERTIES RUN_SERIAL TRUE)
    else()
        _mumps_register_test(${tgt})
    endif()
endfunction()

function(_mumps_attach_link_group tgt variant want_c)
    if(variant STREQUAL "_seq")
        # ``_seq`` variant: libmpiseq stands in for ALL of MPI + BLACS +
        # PBLAS + ScaLAPACK, so the link list is mumps + lapack + blas
        # + multifloats + the test harness. cmake's $<LINK_GROUP:> can't
        # be used because the migrated archives have PUBLIC chains to
        # scalapack / blacs / pblas / pbblas / ptzblas (which we're
        # deliberately replacing with libmpiseq) — naming a subset
        # triggers cmake's cycle detector. Likewise, plain
        # target_link_libraries on these would drag the PUBLIC chain
        # back in (and Intel MPI along the multifloats_mpi PRIVATE chain
        # — STATIC libs propagate PRIVATE deps to consumers). Bypass
        # cmake's link graph entirely: list every required archive
        # inside a single ``--start-group`` block via $<TARGET_FILE:>,
        # add an explicit add_dependencies so cmake builds them, and
        # pull include directories manually from the targets that
        # carry .mod / header search paths.
        set(_archives
            mumps_test_target          # test harness — needs to scan
            mumps_c_bridge             #   alongside mwmumps
            ${LIB_PAIR_PREFIX}mumps
            ${LIB_PAIR_PREFIX}lapack
            ${LIB_PAIR_PREFIX}blas
            mumps_test_common
            target_conv
            mpiseq                     # MPI / BLACS / scalapack stubs
        )
        # Std non-precision archives carry shared helpers (mumps_abort_,
        # qisnan_ for the kind16 EP NaN check, ...). Append whichever
        # exist in this build.
        foreach(_std mumps_common lapack_common blas_common)
            if(TARGET ${_std})
                list(APPEND _archives ${_std})
            endif()
        endforeach()
        # Std blas / lapack supply precision-independent helpers (lsame_,
        # ilaenv_, xerbla_, ieeeck_, dlamch_, …) that the migrated
        # archives call. They used to resolve from refblas_quad /
        # reflapack_quad before those archives' exports were renamed
        # *_quad_; now they need a non-renamed source on the link line.
        foreach(_std blas lapack)
            if(TARGET ${_std})
                list(APPEND _archives ${_std})
            endif()
        endforeach()
        # Ordering archives (pord/metis/scotch) — the migrated solver + bridge
        # + common members reference their symbols (see _mumps_link_group_archives
        # note); include them in the --start-group block, gated on their targets.
        foreach(_ord pord metis esmumps scotch)
            if(TARGET ${_ord})
                list(APPEND _archives ${_ord})
            endif()
        endforeach()
        if(TARGET ${LIB_PAIR_PREFIX}xblas)
            list(APPEND _archives ${LIB_PAIR_PREFIX}xblas)
        endif()
        # Composite-${LIB_PAIR_PREFIX}blas is an INTERFACE under
        # PARALLEL_BLAS=ON; $<TARGET_FILE:...> can't resolve it.
        # Swap to the underlying _serial. See note above near
        # _mumps_link_group_archives.
        if(PARALLEL_BLAS AND TARGET ${LIB_PAIR_PREFIX}blas_serial)
            list(TRANSFORM _archives REPLACE
                "^${LIB_PAIR_PREFIX}blas$" "${LIB_PAIR_PREFIX}blas_serial")
        endif()
        if(NEEDS_MULTIFLOATS)
            # multifloats_mpi (C++ MPI bridge) included as a raw
            # archive — its calls to MPI_Type_create_struct /
            # MPI_Op_create / etc. resolve against libmpiseq's
            # mpiseq_c_stubs.c. The registration runs but produces
            # zeroed handles. Single-rank in-place reductions
            # short-circuit before MUMPS_COPY sees the handle, so
            # this is harmless for MUMPS calls that use MPI_IN_PLACE;
            # non-in-place multifloats reductions would STOP in
            # libseq's MUMPS_COPY (datatype unknown) — that's a real
            # limitation of the libmpiseq path on multifloats and
            # surfaces at runtime if hit.
            list(APPEND _archives
                multifloats_mpi_f multifloats_mpi
                multifloatsf multifloats)
        endif()
        if(NEEDS_QUAD_MPI)
            # quad_mpi (C MPI bridge) + quad_mpi_f (Fortran handles)
            # included as raw archives — quad_mpi's MPI_Op_create calls
            # resolve against libmpiseq's mpiseq_c_stubs.c. Registration
            # runs but produces zeroed handles; single-rank in-place
            # reductions short-circuit before libseq's MUMPS_COPY sees
            # the handle, so this is harmless for the MPI_IN_PLACE path
            # MUMPS uses at np=1. Mirrors the multifloats branch above.
            list(APPEND _archives quad_mpi_f quad_mpi)
        endif()
        # Under the binary-release acquisition the bare ``multifloats``
        # target is an INTERFACE wrapper over the imported
        # ``multifloats::multifloats`` core (cmake/CMakeLists.txt) —
        # $<TARGET_FILE:> can't resolve an INTERFACE library, and
        # add_dependencies() rejects imported targets. Swap in the
        # imported core for the link line and drop it from the
        # dependency list (imported archives are always up to date).
        set(_dep_archives ${_archives})
        if("multifloats" IN_LIST _archives AND TARGET multifloats)
            get_target_property(_mf_type multifloats TYPE)
            if(_mf_type STREQUAL "INTERFACE_LIBRARY")
                list(TRANSFORM _archives REPLACE "^multifloats$"
                     "multifloats::multifloats")
                list(REMOVE_ITEM _dep_archives multifloats)
            endif()
        endif()
        add_dependencies(${tgt} ${_dep_archives})
        # target_link_libraries (vs. target_link_options) places these
        # AFTER the .o files in the link line — required so ld can
        # resolve test_*.o symbols against the group.
        target_link_libraries(${tgt} PRIVATE "-Wl,--start-group")
        foreach(_arc IN LISTS _archives)
            if(_arc STREQUAL "${LIB_PAIR_PREFIX}blas" AND PARALLEL_BLAS
                    AND TARGET ${LIB_PAIR_PREFIX}blas_serial)
                # Composite is INTERFACE → $<TARGET_FILE:> can't
                # resolve. Inline the two underlying archives directly.
                # Overlay first so its egemm_ wins over the serial
                # archive's egemm_ when ld scans the group.
                if(TARGET ${LIB_PAIR_PREFIX}blas_parallel)
                    target_link_libraries(${tgt} PRIVATE
                        $<TARGET_FILE:${LIB_PAIR_PREFIX}blas_parallel>)
                endif()
                target_link_libraries(${tgt} PRIVATE
                    $<TARGET_FILE:${LIB_PAIR_PREFIX}blas_serial>)
            else()
                target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:${_arc}>)
            endif()
        endforeach()
        target_link_libraries(${tgt} PRIVATE "-Wl,--end-group")
        # mumps_test_target / mumps_c_bridge / mumps_test_common /
        # multifloats_mpi_f / target_conv / mpiseq ship .mod files and
        # bridge headers we need at compile time. Pull their public
        # include dirs manually since we're not using
        # target_link_libraries for them.
        foreach(_arc IN LISTS _archives)
            target_include_directories(${tgt} PRIVATE
                $<TARGET_PROPERTY:${_arc},INTERFACE_INCLUDE_DIRECTORIES>)
        endforeach()
        # reflapack_quad / refblas_quad are reference archives the
        # solve helpers reach for. They're optional (built only when
        # tests/lapack and tests/blas are configured).
        if(TARGET reflapack_quad)
            target_link_libraries(${tgt} PRIVATE reflapack_quad)
        endif()
    else()
        set(_archives ${_mumps_link_group_archives})
        target_link_libraries(${tgt} PRIVATE
            "$<LINK_GROUP:RESCAN,${_archives}>"
        )
        # Ordering archives (pord/metis/scotch) appended AFTER the group as raw
        # files. They are leaf libraries — nothing in them calls back into the
        # mumps/scalapack group — so trailing placement resolves the ordering
        # symbols the group members left undefined (ld carries the unresolved
        # set forward). Raw $<TARGET_FILE:> + add_dependencies keeps them out of
        # the target link graph, so the eymumps->mumps_common->pord PUBLIC
        # chain does NOT form a group<->common cycle (which naming them as
        # LINK_GROUP targets would). Their own --start-group covers the
        # esmumps->scotch intra-Scotch cross-references.
        if(_mumps_ordering_archives)
            add_dependencies(${tgt} ${_mumps_ordering_archives})
            target_link_libraries(${tgt} PRIVATE "-Wl,--start-group")
            foreach(_ord IN LISTS _mumps_ordering_archives)
                target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:${_ord}>)
            endforeach()
            target_link_libraries(${tgt} PRIVATE "-Wl,--end-group")
        endif()
    endif()
endfunction()

function(add_mumps_fortran_test test_src)
    get_filename_component(name ${test_src} NAME_WE)
    foreach(_vkey IN LISTS _mumps_test_variants)
        if(_vkey STREQUAL "plain")
            set(variant "")
        else()
            set(variant "${_vkey}")
        endif()
        # PT-Scotch (ICNTL(28)=2/ICNTL(29)=1) is distributed, real-MPI-only
        # analysis: the libptscotch increment is absent from the libseq
        # column, so the *_ptscotch_analysis drivers cannot link (let alone
        # run) under the ``_seq`` (libmpiseq) variant. Build them plain-only.
        if(variant STREQUAL "_seq" AND name MATCHES "ptscotch")
            continue()
        endif()
        # The *_noinit drivers deliberately skip MPI_Init and the per-target
        # MPI-handle registration to prove a plain sequential (libmpiseq)
        # consumer solves with neither. That contract is only meaningful under
        # the ``_seq`` variant; the real-MPI variant legitimately requires
        # MPI_Init, so a no-init driver there would be a false negative. Build
        # them seq-only. (When libmpiseq isn't built ``_seq`` is absent from
        # _mumps_test_variants, so the driver is simply skipped everywhere.)
        if(NOT variant STREQUAL "_seq" AND name MATCHES "noinit")
            continue()
        endif()
        set(tgt ${name}${variant})
        add_executable(${tgt} ${test_src})
        _mumps_attach_link_group(${tgt} "${variant}" FALSE)
        _mumps_apply_link_variant(${tgt} "${variant}" FALSE)
        # Use the native objects in the LTO-fat archives (built with
        # -flto=auto -ffat-lto-objects upstream). Cross-archive LTO
        # partitioning produces spurious "undefined reference to
        # xmumps_mtransi_" / "pxtrsm_" etc. when linking via the LTO
        # plugin, even though those symbols ARE in the archives. Native
        # link sidesteps this.
        target_link_options(${tgt} PRIVATE -fno-lto)
        _mumps_register_test_variant(${tgt} "${variant}")
    endforeach()
endfunction()

function(add_mumps_c_test test_src)
    get_filename_component(name ${test_src} NAME_WE)
    string(TOUPPER "${TARGET_NAME}" _ttn_upper)
    foreach(_vkey IN LISTS _mumps_test_variants)
        if(_vkey STREQUAL "plain")
            set(variant "")
        else()
            set(variant "${_vkey}")
        endif()
        # Seq-only: the *_noinit drivers prove a plain libmpiseq consumer
        # solves without MPI_Init or the per-target MPI-handle registration —
        # a contract that only holds under the ``_seq`` variant. Mirrors the
        # Fortran guard in add_mumps_fortran_test.
        if(NOT variant STREQUAL "_seq" AND name MATCHES "noinit")
            continue()
        endif()
        set(tgt ${name}${variant})
        add_executable(${tgt} ${test_src})
        # Same Fortran-linker / no-LTO setup as the Fortran tests — the
        # bridge archive contains gfortran-LTO bitcode pulled in via
        # libqmumps's LTO partitioning, and ``cc`` doesn't have the
        # gfortran LTO plugin loaded.
        set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
        _mumps_attach_link_group(${tgt} "${variant}" TRUE)
        _mumps_apply_link_variant(${tgt} "${variant}" TRUE)
        target_link_options(${tgt} PRIVATE -fno-lto)
        target_include_directories(${tgt} PRIVATE
            ${_BRIDGE_INC_TARGET} ${_BRIDGE_INC_SHARED} ${_MUMPS_C_INC}
        )
        # Wire the same ``TARGET_*`` macros mumps_c_bridge uses so each
        # test source can ``#include TARGET_REAL_HEADER`` and dispatch
        # to the per-target bridge entry points without hard-coding
        # kind16.
        target_compile_definitions(${tgt} PRIVATE
            "TARGET_REAL_HEADER=\"${LIB_PREFIX}mumps_c.h\""
            "TARGET_COMPLEX_HEADER=\"${LIB_PREFIX_COMPLEX}mumps_c.h\""
            TARGET_REAL_MUMPS_C=${LIB_PREFIX}mumps_c
            TARGET_COMPLEX_MUMPS_C=${LIB_PREFIX_COMPLEX}mumps_c
            TARGET_REAL_STRUC_C=${PFX}MUMPS_STRUC_C
            TARGET_COMPLEX_STRUC_C=${PFX_C}MUMPS_STRUC_C
            "TEST_TARGET_NAME=\"${TARGET_NAME}\""
            "TEST_TARGET_${_ttn_upper}"
        )
        _mumps_register_test_variant(${tgt} "${variant}")
    endforeach()
endfunction()

file(GLOB _ftests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_*.f90)
foreach(_t ${_ftests})
    add_mumps_fortran_test(${_t})
endforeach()

# Per-target fypp Fortran tests (parity tests, etc.). Each .fypp file
# is processed with ``-Dtarget=${TARGET_NAME}`` so the body can switch
# bind(C) interfaces and bridge buffers per precision.
file(GLOB _ffypp_tests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fortran/test_*.fypp)
foreach(_t ${_ffypp_tests})
    get_filename_component(_stem ${_t} NAME_WE)
    set(_out ${CMAKE_CURRENT_BINARY_DIR}/${_stem}.f90)
    # fypp -D evaluates the RHS as Python, so wrap the target name in
    # quotes so it parses as a string literal rather than an identifier.
    add_custom_command(
        OUTPUT ${_out}
        COMMAND ${FYPP} -Dtarget=\"${TARGET_NAME}\" ${_t} ${_out}
        DEPENDS ${_t}
        COMMENT "fypp ${_t} -> ${_out}"
        VERBATIM
    )
    add_mumps_fortran_test(${_out})
endforeach()

# C-side tests now portable across all targets — the multifloats
# branch in test_real_compat.h widens at the bridge boundary.
#
# test_mumps_c_coexist.c is NOT an add_mumps_c_test: it is a link-only
# probe against the SHIPPED bridge archives (not the test-tree
# mumps_c_bridge), wired separately just below. Exclude it from the glob.
file(GLOB _ctests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/c/test_*.c)
list(FILTER _ctests EXCLUDE REGEX "/test_mumps_c_coexist\\.c$")
foreach(_t ${_ctests})
    add_mumps_c_test(${_t})
endforeach()

# --- C bridge coexistence link checks --------------------------------
# A release emits FOUR ?mumps_c entry points: the migrated real/complex
# pair folded into ${LIB_PAIR_PREFIX}mumps (e/y, q/x or m/w — always the
# upstream d/z lineages) and the genuine double pair folded into dzmumps.
# All share the arithmetic-independent C runtime, now folded into
# libmumps_common.
#
# test_mumps_c_coexist.c names all four and is force-linked (whole
# archive) against each COHERENT family; each MUST link with no
# "multiple definition" of a shared runtime symbol. The two families
# are deliberately never linked together: the migrated real bridge is
# compiled at MUMPS_ARITH_d and so emits the same dmumps_* Fortran
# callbacks as the genuine d bridge (likewise z), and their Fortran
# commons collide — one family per program, exactly as upstream ships
# one libmumps_common per program. Unresolved Fortran/MPI/cross-family
# symbols are ignored: this is a pure duplicate-symbol link check, not a
# run. The tests drive the build of an EXCLUDE_FROM_ALL target, so a
# duplicate-symbol regression fails only this test, loudly.
set(_coexist_src ${CMAKE_CURRENT_SOURCE_DIR}/c/test_mumps_c_coexist.c)
set(_coexist_defs
    TARGET_REAL_MUMPS_C=${LIB_PREFIX}mumps_c
    TARGET_COMPLEX_MUMPS_C=${LIB_PREFIX_COMPLEX}mumps_c)

# genuine double pair — the capability the shared runtime unlocks
# (standard upstream "use dmumps_c and zmumps_c in one program").
if(TARGET dzmumps AND TARGET mumps_common)
    add_executable(test_mumps_c_coexist_genuine EXCLUDE_FROM_ALL ${_coexist_src})
    target_compile_definitions(test_mumps_c_coexist_genuine PRIVATE ${_coexist_defs})
    target_link_options(test_mumps_c_coexist_genuine PRIVATE
        -Wl,--whole-archive
        $<TARGET_FILE:dzmumps> $<TARGET_FILE:mumps_common>
        -Wl,--no-whole-archive
        -Wl,--unresolved-symbols=ignore-all
        -fno-lto)
    add_dependencies(test_mumps_c_coexist_genuine dzmumps mumps_common)
    add_test(NAME mumps_c_coexist_genuine
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                --target test_mumps_c_coexist_genuine --config $<CONFIG>)
endif()

# migrated pair — already shipped; guards the ${LIB_PAIR_PREFIX}mumps bridge
# (both e/y-lineage objects, now folded in) against a future
# runtime-duplication regression.
if(TARGET ${LIB_PAIR_PREFIX}mumps AND TARGET mumps_common)
    add_executable(test_mumps_c_coexist_migrated EXCLUDE_FROM_ALL ${_coexist_src})
    target_compile_definitions(test_mumps_c_coexist_migrated PRIVATE ${_coexist_defs})
    target_link_options(test_mumps_c_coexist_migrated PRIVATE
        -Wl,--whole-archive
        $<TARGET_FILE:${LIB_PAIR_PREFIX}mumps> $<TARGET_FILE:mumps_common>
        -Wl,--no-whole-archive
        -Wl,--unresolved-symbols=ignore-all
        -fno-lto)
    add_dependencies(test_mumps_c_coexist_migrated ${LIB_PAIR_PREFIX}mumps mumps_common)
    add_test(NAME mumps_c_coexist_migrated
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}
                --target test_mumps_c_coexist_migrated --config $<CONFIG>)
endif()

# --- Genuine double-precision solve test -----------------------------
# The d/z-named tests under c/ are remapped (LIB_PREFIX / LIB_PREFIX_COMPLEX)
# onto the MIGRATED bridge (q/x, e/y or m/w) and so NEVER exercise the
# genuine double dmumps_c / zmumps_c — coexist_genuine above is only a
# link probe (--unresolved-symbols=ignore-all, never run). This test
# fills that gap: one executable runs a real solve through genuine
# dmumps_c and a complex solve through genuine zmumps_c, each to
# plain-double accuracy (~1e-16), host-guarded so it is correct under
# mpiexec -n>1 as well as MPISEQ. Kept under c/genuine/ so the
# c/test_*.c glob (compiled against the migrated bridge) skips it.
#
# Genuine and migrated collide on the unprefixed Fortran commons, so the
# genuine link group is a DISJOINT archive set from the migrated one:
# dzmumps (with its folded-in C bridge) + mumps_common (Fortran common +
# the folded-in shared C runtime) + the MPI-matched in-tree parallel stack
# (cmake/CMakeLists.txt's _genuine_parallel_stack: bundled
# ScaLAPACK/PBLAS/BLACS/LAPACK/BLAS built against the chosen MPI). The
# circular static refs among those archives are resolved by the same
# LINK_GROUP RESCAN the migrated tests use. The _seq variant swaps the
# whole parallel stack for libmpiseq, which natively supplies the
# DOUBLE / PD* BLACS + ScaLAPACK stubs the genuine solvers call.
set(_genuine_solve_src ${CMAKE_CURRENT_SOURCE_DIR}/c/genuine/test_genuine_dz_solve.c)
if(TARGET dzmumps
        AND TARGET mumps_common AND EXISTS ${_genuine_solve_src})
    # "" variant: genuine solvers (bridge folded in) + common (with the
    # folded-in C runtime) + full parallel stack.
    set(_genuine_group
        dzmumps mumps_common)
    foreach(_p scalapack scalapack_common pblas pblas_common pbblas
               ptzblas ptzblas_common blacs blacs_common
               lapack lapack_common blas)
        if(TARGET ${_p})
            list(APPEND _genuine_group ${_p})
        endif()
    endforeach()
    # Ordering archives (ICNTL(7)=4 pord, =5 metis, =3 scotch). mumps_common
    # PUBLIC-links these, but its transitive interface is NOT pulled *into* the
    # LINK_GROUP RESCAN below, so the ordering symbols the group members reference
    # (mumps_pord.c's SPACE_ordering/firstPostorder, mumps_metis.c's
    # metis_mumps_*, mumps_scotch.c's SCOTCH_*_mumps, and the Fortran
    # ana_orderings_wrappers_m / d/z ini_driver callers) would be undefined at
    # link. Name them explicitly inside the group, gated on the target so the
    # tree still builds when a given ordering wasn't staged.
    foreach(_p pord metis esmumps scotch)
        if(TARGET ${_p})
            list(APPEND _genuine_group ${_p})
        endif()
    endforeach()
    # PT-Scotch increment for the genuine (MKL d/z) stack: dzmumps is built
    # with -Dptscotch too (cmake/CMakeLists.txt), so it references the
    # distributed SCOTCH_dgraph*_mumps / SCOTCHFDGRAPH*_mumps symbols. MPI-only
    # → plain variant only; gated on the target so libseq stays sequential.
    if(TARGET ptscotch)
        list(APPEND _genuine_group ptscotch)
    endif()
    # _seq variant: libmpiseq replaces ScaLAPACK/PBLAS/BLACS entirely;
    # keep the std double LAPACK + BLAS the genuine solvers call.
    set(_genuine_group_seq
        dzmumps mumps_common)
    foreach(_p lapack blas)
        if(TARGET ${_p})
            list(APPEND _genuine_group_seq ${_p})
        endif()
    endforeach()
    # Same ordering archives for the sequential stack (mumps_common is present
    # in both variants and references them identically).
    foreach(_p pord metis esmumps scotch)
        if(TARGET ${_p})
            list(APPEND _genuine_group_seq ${_p})
        endif()
    endforeach()

    foreach(_vkey IN LISTS _mumps_test_variants)
        if(_vkey STREQUAL "plain")
            set(variant "")
        else()
            set(variant "${_vkey}")
        endif()
        set(tgt test_genuine_dz_solve${variant})
        add_executable(${tgt} ${_genuine_solve_src})
        # Fortran linker + no-LTO: the genuine bridges pull gfortran-LTO
        # bitcode from libdmumps/libzmumps, exactly like the migrated
        # C tests (see add_mumps_c_test).
        set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
        target_link_options(${tgt} PRIVATE -fno-lto)
        # Upstream dmumps_c.h / zmumps_c.h / mumps_c_types.h.
        target_include_directories(${tgt} PRIVATE ${_MUMPS_C_INC})
        if(variant STREQUAL "_seq")
            add_dependencies(${tgt} ${_genuine_group_seq} mpiseq)
            target_link_libraries(${tgt} PRIVATE "-Wl,--start-group")
            foreach(_arc IN LISTS _genuine_group_seq)
                target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:${_arc}>)
            endforeach()
            target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:mpiseq>)
            target_link_libraries(${tgt} PRIVATE "-Wl,--end-group")
            foreach(_arc IN LISTS _genuine_group_seq mpiseq)
                target_include_directories(${tgt} PRIVATE
                    $<TARGET_PROPERTY:${_arc},INTERFACE_INCLUDE_DIRECTORIES>)
            endforeach()
        else()
            target_link_libraries(${tgt} PRIVATE
                "$<LINK_GROUP:RESCAN,${_genuine_group}>")
        endif()
        _mumps_apply_link_variant(${tgt} "${variant}" TRUE)
        _mumps_register_test_variant(${tgt} "${variant}")
    endforeach()

    # ── genuine Fortran-API solve test ───────────────────────────────
    # Companion to the C-bridge test above: calls the genuine DMUMPS /
    # ZMUMPS Fortran entry points directly on the DMUMPS_STRUC /
    # ZMUMPS_STRUC derived types (upstream's canonical usage). The C
    # bridge only reaches libdmumps/libzmumps through the dmumps_f77_
    # shim, so the Fortran calling convention and derived-type layout are
    # covered ONLY here. Reuses the same disjoint genuine link groups
    # (the solver needs mumps_common for its C ordering/IO symbols);
    # want_c is FALSE — a pure Fortran executable, no C bridge object.
    set(_genuine_solve_f_src
        ${CMAKE_CURRENT_SOURCE_DIR}/fortran/genuine/test_genuine_dz_solve_f.F90)
    if(EXISTS ${_genuine_solve_f_src})
        foreach(_vkey IN LISTS _mumps_test_variants)
            if(_vkey STREQUAL "plain")
                set(variant "")
            else()
                set(variant "${_vkey}")
            endif()
            set(tgt test_genuine_dz_solve_f${variant})
            add_executable(${tgt} ${_genuine_solve_f_src})
            target_link_options(${tgt} PRIVATE -fno-lto)
            # dmumps_struc.h / zmumps_struc.h live in the upstream include.
            target_include_directories(${tgt} PRIVATE ${_MUMPS_C_INC})
            if(variant STREQUAL "_seq")
                add_dependencies(${tgt} ${_genuine_group_seq} mpiseq)
                target_link_libraries(${tgt} PRIVATE "-Wl,--start-group")
                foreach(_arc IN LISTS _genuine_group_seq)
                    target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:${_arc}>)
                endforeach()
                target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:mpiseq>)
                target_link_libraries(${tgt} PRIVATE "-Wl,--end-group")
                foreach(_arc IN LISTS _genuine_group_seq mpiseq)
                    target_include_directories(${tgt} PRIVATE
                        $<TARGET_PROPERTY:${_arc},INTERFACE_INCLUDE_DIRECTORIES>)
                endforeach()
            else()
                target_link_libraries(${tgt} PRIVATE
                    "$<LINK_GROUP:RESCAN,${_genuine_group}>")
            endif()
            _mumps_apply_link_variant(${tgt} "${variant}" FALSE)
            _mumps_register_test_variant(${tgt} "${variant}")
        endforeach()
    endif()

    # ── genuine ALL-precision C-bridge solve test (s + c + d + z) ─────
    # The user-requested "link all precision libs into one executable
    # which tests them all": libsmumps / libcmumps / libdmumps /
    # libzmumps coexist in a single program because they share the one
    # unpromoted mumps_common (what retiring dzmumps_common unlocked) and
    # export distinct {s,c,d,z}mumps_c entry points. Extends the genuine
    # link group with the single-precision solver + bridge archives; the
    # migrated e/y stack stays out (its bridge clashes on dmumps_* C
    # helpers — see the file header). Guarded on all four genuine bridges
    # so it only appears in the genuine-enabled build.
    set(_genuine_scdz_src
        ${CMAKE_CURRENT_SOURCE_DIR}/c/genuine/test_genuine_scdz_solve.c)
    if(TARGET scmumps
            AND EXISTS ${_genuine_scdz_src})
        # scmumps solver (bridge folded in) prepended to the genuine d/z
        # groups; the shared common + parallel stack are already listed there.
        set(_genuine_scdz_group     scmumps ${_genuine_group})
        set(_genuine_scdz_group_seq scmumps ${_genuine_group_seq})
        foreach(_vkey IN LISTS _mumps_test_variants)
            if(_vkey STREQUAL "plain")
                set(variant "")
            else()
                set(variant "${_vkey}")
            endif()
            set(tgt test_genuine_scdz_solve${variant})
            add_executable(${tgt} ${_genuine_scdz_src})
            set_target_properties(${tgt} PROPERTIES LINKER_LANGUAGE Fortran)
            target_link_options(${tgt} PRIVATE -fno-lto)
            target_include_directories(${tgt} PRIVATE ${_MUMPS_C_INC})
            if(variant STREQUAL "_seq")
                add_dependencies(${tgt} ${_genuine_scdz_group_seq} mpiseq)
                target_link_libraries(${tgt} PRIVATE "-Wl,--start-group")
                foreach(_arc IN LISTS _genuine_scdz_group_seq)
                    target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:${_arc}>)
                endforeach()
                target_link_libraries(${tgt} PRIVATE $<TARGET_FILE:mpiseq>)
                target_link_libraries(${tgt} PRIVATE "-Wl,--end-group")
                foreach(_arc IN LISTS _genuine_scdz_group_seq mpiseq)
                    target_include_directories(${tgt} PRIVATE
                        $<TARGET_PROPERTY:${_arc},INTERFACE_INCLUDE_DIRECTORIES>)
                endforeach()
            else()
                target_link_libraries(${tgt} PRIVATE
                    "$<LINK_GROUP:RESCAN,${_genuine_scdz_group}>")
            endif()
            _mumps_apply_link_variant(${tgt} "${variant}" TRUE)
            _mumps_register_test_variant(${tgt} "${variant}")
        endforeach()
    endif()
endif()
