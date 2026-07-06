# eplinalg — MUMPS: orderings, typed solvers, C bridge, genuine dz/sc
#
# include()'d from the staged top-level CMakeLists.txt (same directory
# scope — variables and targets defined here behave exactly as if the
# content were inline). Staged flat next to CMakeLists.txt by
# ``migrator stage`` (see the cmake-glue copy lists in
# ``src/migrator/staging.py``).

# 9. MUMPS (depends on LAPACK, BLAS, ScaLAPACK; needs MPI Fortran).
#
# MUMPS is a sparse direct solver (Multifrontal Massively Parallel). The
# migrated archive ``${LIB_PAIR_PREFIX}mumps`` carries the Q/X clones of
# DMUMPS / ZMUMPS plus their auxiliary subroutines.
#
# (The genuine-precision solvers are handled separately at the end of
# this section: libdzmumps / libscmumps rebuild the upstream s/c/d/z
# arithmetics against the in-tree parallel stack.)
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
            # Ship the archive as libpord_mumps.a (not libpord.a): the
            # CMake target keeps its short name ``pord`` for the internal
            # link graph, but the installed filename carries the _mumps
            # tag so it cannot collide on disk with a system/other-dev
            # ordering archive of the same base name (mirrors the _mumps
            # symbol namespacing).
            set_target_properties(pord PROPERTIES OUTPUT_NAME pord_mumps)
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
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON
                # libmetis_mumps.a — filename-namespaced like the _mumps
                # symbol rename, so it never clashes with a system METIS.
                OUTPUT_NAME metis_mumps)
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
    # gains no bison/flex dependency. Building ``scotch`` (+ ``esmumps``,
    # the MUMPS adapter) here and defining ``-Dscotch`` on the
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

            # Ship these two as lib{scotch,esmumps}_mumps.a: the _mumps
            # filename tag mirrors the SCOTCH_NAME_SUFFIX=_mumps symbol
            # namespacing, so this private Scotch copy cannot clash on
            # disk with a system Scotch another developer might add. The
            # CMake target names stay short for the internal link graph.
            #
            # Upstream ships the error handler (library_error*.c) as a
            # separate libscotcherr.a so applications can substitute
            # their own handler archive on the link line. We fold it
            # into libscotch itself: ld extracts archive members on
            # demand, so an application-supplied handler *object* still
            # overrides the archive member, and nothing here swaps
            # handlers — one archive fewer to ship and order.
            add_library(scotch STATIC ${_scotch_lib_c} ${_scotch_err_c})
            target_include_directories(scotch PUBLIC
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>
                $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
            target_compile_definitions(scotch PRIVATE ${_scotch_defs})
            set_target_properties(scotch PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON
                OUTPUT_NAME scotch_mumps)

            add_library(esmumps STATIC ${_scotch_esm_c})
            target_include_directories(esmumps PUBLIC
                $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
            target_include_directories(esmumps PRIVATE
                $<BUILD_INTERFACE:${_mumps_scotch_esmsrc}>
                $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
            target_compile_definitions(esmumps PRIVATE ${_scotch_defs})
            set_target_properties(esmumps PROPERTIES
                C_STANDARD 99 POSITION_INDEPENDENT_CODE ON
                OUTPUT_NAME esmumps_mumps)
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
                target_compile_options(scotch  PRIVATE ${_scotch_cflags})
                target_compile_options(esmumps PRIVATE ${_scotch_cflags})
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

    # ── PT-Scotch (distributed parallel analysis) ───────────────────
    # The 131-file distributed increment (dgraph*/bdgraph*/vdgraph*/
    # hdgraph*/kdgraph*/dorder*/dmapping*/library_dgraph*/library_context_
    # dgraph*, from SCOTCH_LIBPTSCOTCH_SOURCES). Compiled WITH
    # -DSCOTCH_PTSCOTCH (which chains module.h → COMMON_MPI → <mpi.h>) and
    # the SAME -DSCOTCH_NAME_SUFFIX=_mumps rename as sequential Scotch, so
    # the archive is a disjoint increment that references the sequential
    # symbols by their _mumps names and resolves them from libscotch_mumps
    # at link (upstream's -lptscotch -lscotch model). It is genuinely
    # MPI-tied, so it is gated on MPI and installs -${MPI_TAG}-tagged
    # (_MPI_DEPENDENT_LIBS below); the libseq column keeps it off. Turning
    # this on plus -Dptscotch on the MUMPS runtime activates the
    # ICNTL(28)=2 + ICNTL(29)=1 distributed-analysis path, complementing
    # the sequential ICNTL(7)=3 Scotch that libscotch already serves.
    #
    # ptscotch.h / ptscotchf.h are pre-generated & vendored in the include
    # dir (like scotch.h/scotchf.h): dummysizes computes the opaque-struct
    # dims, and those dims are provably MPI-ABI-independent — Dgraph embeds
    # an MPI_Comm by value (8 bytes on OpenMPI, 4 on MPICH/Intel-MPI), but
    # the struct is 8-byte aligned (it holds pointers), so sizeof is always
    # a multiple of 8 and the ≤4-byte delta rounds back up. Verified: the
    # OpenMPI- and Intel-MPI-generated headers are byte-identical.
    # PT-Scotch is a distributed (real-MPI) ordering; the libmpiseq (seq)
    # release links no real MPI, so it is forced off there — matching the
    # "libseq column keeps PT-Scotch off" rule and keeping mumps_common's
    # -Dptscotch out of a link line that has no ptscotch archive.
    set(MUMPS_HAVE_PTSCOTCH FALSE)
    if(MUMPS_HAVE_SCOTCH AND MPI_C_FOUND AND MPI_Fortran_FOUND
            AND SCOTCH_LIBPTSCOTCH_SOURCES AND NOT MUMPS_LIBSEQ_RELEASE)
        set(_ptscotch_lib_c "")
        foreach(_f ${SCOTCH_LIBPTSCOTCH_SOURCES})
            list(APPEND _ptscotch_lib_c ${_mumps_scotch_libsrc}/${_f})
        endforeach()
        # Sequential defines + the distributed flag. The suffix MUST stay
        # _mumps (not _mumps_pt) — see increment rationale above.
        set(_ptscotch_defs ${_scotch_defs} SCOTCH_PTSCOTCH)

        # Filename carries the MPI tag (unlike the shared, untagged
        # sequential libscotch_mumps.a): OpenMPI and Intel-MPI copies are
        # different binaries and must coexist in one lib/ — same discipline
        # as libdmumps-${MPI_TAG}.a. Falls back to plain _mumps if MPI_TAG
        # is empty (unknown vendor). Baked in here rather than via
        # fortran_install_library, which would drop the _mumps stem and add
        # the compiler tag.
        set(_ptscotch_out ptscotch_mumps)
        if(MPI_LIB_TAG)
            set(_ptscotch_out ptscotch_mumps-${MPI_LIB_TAG})
        endif()

        # The error handler (library_error*.c) is rebuilt under
        # -DSCOTCH_PTSCOTCH — a distinct, MPI-aware (rank-printing)
        # object from the one folded into sequential libscotch — and
        # folded into libptscotch the same way. On a pt link line
        # (ptscotch before scotch) the distributed increment's error
        # references extract THIS member first, so the seq copy inside
        # libscotch is never pulled — the MPI-aware handler wins with
        # no duplicate definitions, mirroring upstream's -lptscotcherr
        # -lscotch -lscotcherr order with two fewer archives.
        add_library(ptscotch STATIC ${_ptscotch_lib_c} ${_scotch_err_c})
        target_include_directories(ptscotch PUBLIC
            $<BUILD_INTERFACE:${_mumps_scotch_inc}>
            $<BUILD_INTERFACE:${_mumps_scotch_libsrc}>)
        target_compile_definitions(ptscotch PRIVATE ${_ptscotch_defs})
        # The distributed increment references the sequential SCOTCH_*_mumps
        # symbols but does not define them. Declaring the dependency PUBLIC
        # makes CMake emit ptscotch BEFORE scotch on every consumer's static
        # link line (GNU ld is single-pass), so those forward references
        # resolve — mirroring upstream's -lptscotch -lscotch order. Every
        # consumer then only has to link ptscotch; scotch follows
        # transitively in the correct order and dedups against any explicit
        # scotch entry.
        target_link_libraries(ptscotch PUBLIC scotch MPI::MPI_C)
        set_target_properties(ptscotch PROPERTIES
            C_STANDARD 99 POSITION_INDEPENDENT_CODE ON
            OUTPUT_NAME ${_ptscotch_out})

        # Same warning-quiet + based-array fortify mitigation as sequential
        # Scotch (the distributed sources use the identical idiom).
        if(NOT MSVC)
            target_compile_options(ptscotch PRIVATE ${_scotch_cflags})
        endif()

        # Distributed Fortran entry points MUMPS calls under -Dptscotch
        # ([dscz]ana_aux_par.F, ana_orderings_wrappers_m.F). Same
        # suffix-before-trailing-underscore convention as _scotch_f_renames.
        # SCOTCHFSTRATINIT/EXIT are shared base symbols already remapped by
        # _scotch_f_renames (applied together, since scotch is always on
        # when ptscotch is), and DGRAPH init is done via the C shim
        # (SCOTCH_dgraphInit), so neither appears here. The CONTEXT* set is
        # included because dana_aux_par.F names them (some behind
        # MUMPS_SCOTCHIMPORTOMPTHREADS, some — CONTEXTBINDDGRAPH — under a
        # plain #if defined(ptscotch)); a rename for an unreferenced symbol
        # is inert.
        set(_scotch_pt_f_renames
            SCOTCHFDGRAPHBUILD=SCOTCHFDGRAPHBUILD_mumps
            SCOTCHFDGRAPHEXIT=SCOTCHFDGRAPHEXIT_mumps
            SCOTCHFDGRAPHORDERINIT=SCOTCHFDGRAPHORDERINIT_mumps
            SCOTCHFDGRAPHORDEREXIT=SCOTCHFDGRAPHORDEREXIT_mumps
            SCOTCHFDGRAPHORDERCOMPUTE=SCOTCHFDGRAPHORDERCOMPUTE_mumps
            SCOTCHFDGRAPHORDERGATHER=SCOTCHFDGRAPHORDERGATHER_mumps
            SCOTCHFDGRAPHCORDERINIT=SCOTCHFDGRAPHCORDERINIT_mumps
            SCOTCHFDGRAPHCORDEREXIT=SCOTCHFDGRAPHCORDEREXIT_mumps
            SCOTCHFSTRATDGRAPHORDER=SCOTCHFSTRATDGRAPHORDER_mumps
            SCOTCHFCONTEXTINIT=SCOTCHFCONTEXTINIT_mumps
            SCOTCHFCONTEXTEXIT=SCOTCHFCONTEXTEXIT_mumps
            SCOTCHFCONTEXTRANDOMCLONE=SCOTCHFCONTEXTRANDOMCLONE_mumps
            SCOTCHFCONTEXTBINDDGRAPH=SCOTCHFCONTEXTBINDDGRAPH_mumps
            SCOTCHFCONTEXTTHREADIMPORT1=SCOTCHFCONTEXTTHREADIMPORT1_mumps
            SCOTCHFCONTEXTTHREADIMPORT2=SCOTCHFCONTEXTTHREADIMPORT2_mumps)

        set(MUMPS_HAVE_PTSCOTCH TRUE)
        message(STATUS "MUMPS: PT-Scotch parallel analysis enabled "
                       "(ICNTL(28)=2, ICNTL(29)=1)")
    endif()
    if(NOT MUMPS_HAVE_PTSCOTCH)
        if(MUMPS_LIBSEQ_RELEASE)
            message(STATUS "MUMPS: PT-Scotch off — libmpiseq (seq) release "
                           "links no real MPI; ICNTL(28)=2 inactive")
        else()
            message(STATUS "MUMPS: PT-Scotch parallel analysis unavailable "
                           "(needs Scotch + MPI C/Fortran) — ICNTL(28)=2 inactive")
        endif()
    endif()

    # _add_mumps_ordering_defines(<target> [FORTRAN_ONLY]): put the
    # ordering ``-D`` defines (-Dpord/-Dmetis/-Dscotch/-Dptscotch, each
    # gated on its MUMPS_HAVE_* flag) on a MUMPS Fortran target. The
    # analysis Fortran (ana_set_ordering.F, ana_orderings_wrappers_m.F,
    # *ana_aux*.F, tools_common.F, mumps_print_defined.F) compiles its
    # ordering call paths only under these defines; without them the
    # branches drop to the "ordering unavailable" error (PT-Scotch:
    # INFO(1)=-38). The referenced entry points resolve from
    # mumps_common's folded C runtime (MUMPS_PORDF / METIS_MUMPS_* /
    # SCOTCH_*_mumps) and the ordering archives it links. Scotch adds
    # the CPP renames binding bare SCOTCHF*/CONTEXT* Fortran calls to
    # the _mumps-suffixed archive symbols, plus the include dir for
    # INCLUDE 'scotchf.h' (ptscotchf.h shares it); PT-Scotch adds the
    # distributed rename set. FORTRAN_ONLY gates the plain defines to
    # Fortran (for targets whose C sources must not see them); the
    # renames are always Fortran-gated. Ordering-archive LINKS stay at
    # the call sites — each target wires them differently.
    function(_add_mumps_ordering_defines target)
        cmake_parse_arguments(_MOD "FORTRAN_ONLY" "" "" ${ARGN})
        foreach(_ord pord metis scotch ptscotch)
            string(TOUPPER "${_ord}" _ord_up)
            if(NOT MUMPS_HAVE_${_ord_up})
                continue()
            endif()
            if(_MOD_FORTRAN_ONLY)
                target_compile_definitions(${target} PRIVATE
                    $<$<COMPILE_LANGUAGE:Fortran>:${_ord}>)
            else()
                target_compile_definitions(${target} PRIVATE ${_ord})
            endif()
            if(_ord STREQUAL "scotch")
                target_compile_definitions(${target} PRIVATE
                    $<$<COMPILE_LANGUAGE:Fortran>:${_scotch_f_renames}>)
                target_include_directories(${target} PRIVATE
                    $<BUILD_INTERFACE:${_mumps_scotch_inc}>)
            elseif(_ord STREQUAL "ptscotch")
                target_compile_definitions(${target} PRIVATE
                    $<$<COMPILE_LANGUAGE:Fortran>:${_scotch_pt_f_renames}>)
            endif()
        endforeach()
    endfunction()

    if(TARGET ${LIB_PAIR_PREFIX}mumps)
        foreach(_dep ${LIB_PAIR_PREFIX}scalapack ${LIB_PAIR_PREFIX}lapack ${LIB_PAIR_PREFIX}blas)
            if(TARGET ${_dep})
                target_link_libraries(${LIB_PAIR_PREFIX}mumps PUBLIC ${_dep})
            endif()
        endforeach()
        if(MPI_Fortran_FOUND)
            # PRIVATE so consumers can substitute libmpiseq at link time
            # without dragging libmpifort.so along the PUBLIC dep chain.
            target_link_libraries(${LIB_PAIR_PREFIX}mumps PRIVATE MPI::MPI_Fortran)
        endif()
        # Multifloats: the migrator rewrites MPI_DOUBLE_PRECISION →
        # MPI_FLOAT64X2 (and reduction op tokens) in MUMPS Fortran sources,
        # then injects ``USE multifloats_mpi_f`` so those names resolve.
        # Link the module library so the .mod is on the include path and
        # the C++ side that actually defines the handles is pulled in.
        if(NEEDS_MULTIFLOATS AND TARGET multifloats_mpi_f)
            target_link_libraries(${LIB_PAIR_PREFIX}mumps PUBLIC multifloats_mpi_f)
        endif()
        # Quad (kind16): same story — the migrator rewrites the reduce-op
        # tokens to MPI_QQ_SUM / MPI_XX_* on the standard MPI_REAL16 /
        # MPI_COMPLEX32 datatypes and injects ``USE quad_mpi_f``. Link the
        # shim so its .mod is on the include path and quad_mpi (which
        # defines the handles) is pulled in transitively.
        if(NEEDS_QUAD_MPI AND TARGET quad_mpi_f)
            target_link_libraries(${LIB_PAIR_PREFIX}mumps PUBLIC quad_mpi_f)
        endif()
        # Ordering defines on the migrated analysis Fortran; the
        # defines reach it because CMAKE_Fortran_PREPROCESS is ON.
        # Link the distributed increment BEFORE the sequential library
        # (upstream -lptscotch -lscotch order; scotch itself is reached
        # via mumps_common's folded C runtime).
        _add_mumps_ordering_defines(${LIB_PAIR_PREFIX}mumps)
        if(MUMPS_HAVE_PTSCOTCH)
            target_link_libraries(${LIB_PAIR_PREFIX}mumps PUBLIC ptscotch)
        endif()
    endif()
    if(TARGET mumps_common AND MPI_Fortran_FOUND)
        target_link_libraries(mumps_common PRIVATE MPI::MPI_Fortran)
    endif()
    # Same defines on the common analysis Fortran (ana_set_ordering.F /
    # ana_orderings_wrappers_m.F / tools_common.F / mumps_print_defined.F
    # gate their ordering branches on them; they apply to the folded C
    # runtime too). The C-side force-includes and ordering-archive links
    # are added in the C-runtime block below.
    if(TARGET mumps_common)
        _add_mumps_ordering_defines(mumps_common)
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
    # Built here (not test-only) so ``find_package(qxmumps)`` delivers a
    # usable C API: the raw upstream ``dmumps_c.h`` alone would be
    # misleading — it expects libdmumps, not libqmumps.
    set(_mumps_c_src ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_src)
    set(_mumps_c_inc ${CMAKE_CURRENT_SOURCE_DIR}/_mumps_upstream_include)
    set(_mumps_c_bridge_inc_target
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/target_${TARGET_NAME}/c/include)
    set(_mumps_c_bridge_inc_shared
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/mumps/c/include)

    # _add_mumps_c_object(<name> <arith> [EXTENDED <prefix>]):
    # compile upstream mumps_c.c once as OBJECT library <name> for
    # arithmetic <arith> (one of s/c/d/z). Without EXTENDED it is a
    # genuine native-width bridge: the upstream include dir comes FIRST
    # so mumps_c_types.h keeps its plain float/double widths. With
    # EXTENDED <prefix> the per-target bridge include dirs come first,
    # mumps_c_types_extended.h is force-included (overriding the double
    # widths with the target's extended types), and the upstream
    # <arith>mumps_* entry/struct names are macro-renamed onto the
    # <prefix>mumps_* symbols so both arithmetics of a migrated family
    # land in one archive without clashing.
    function(_add_mumps_c_object name arith)
        cmake_parse_arguments(_MCO "" "EXTENDED" "" ${ARGN})
        add_library(${name} OBJECT ${_mumps_c_src}/mumps_c.c)
        if(_MCO_EXTENDED)
            string(TOUPPER "${arith}" _arith_up)
            string(TOUPPER "${_MCO_EXTENDED}" _prefix_up)
            target_include_directories(${name} PRIVATE
                ${_mumps_c_bridge_inc_target} ${_mumps_c_bridge_inc_shared}
                ${_mumps_c_inc} ${_mumps_c_src})
            target_compile_definitions(${name} PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_${arith} Add_
                ${arith}mumps_f77_=${_MCO_EXTENDED}mumps_f77_
                ${arith}mumps_set_tmp_ptr_=${_MCO_EXTENDED}mumps_set_tmp_ptr_
                ${arith}mumps_c=${_MCO_EXTENDED}mumps_c
                ${_arith_up}MUMPS_STRUC_C=${_prefix_up}MUMPS_STRUC_C)
            target_compile_options(${name} PRIVATE
                -include mumps_c_types_extended.h)
        else()
            target_include_directories(${name} PRIVATE
                ${_mumps_c_inc} ${_mumps_c_src} ${_mumps_c_bridge_inc_shared})
            target_compile_definitions(${name} PRIVATE
                MUMPS_ARITH=MUMPS_ARITH_${arith} Add_)
        endif()
    endfunction()
    if(TARGET ${LIB_PAIR_PREFIX}mumps
            AND IS_DIRECTORY ${_mumps_c_src}
            AND IS_DIRECTORY ${_mumps_c_inc}
            AND IS_DIRECTORY ${_mumps_c_bridge_inc_target})
        # Per-arithmetic OBJECT libraries — compile ``mumps_c.c`` twice
        # with per-arith macro renames; the ``mumps_c.h`` ↔
        # ``dmumps_c.h``-style override gives us ``qmumps_c`` (real)
        # and ``xmumps_c`` (complex) entry points landing in the same
        # archive.
        _add_mumps_c_object(${LIB_PREFIX}mumps_c_obj d
            EXTENDED ${LIB_PREFIX})
        _add_mumps_c_object(${LIB_PREFIX_COMPLEX}mumps_c_obj z
            EXTENDED ${LIB_PREFIX_COMPLEX})

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
        # duplicate-symbol clash. The C sources keep their required
        # include order + force-include via C-language-gated per-target
        # options below.
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
        # (${LIB_PAIR_PREFIX}mumps PUBLIC-links it via add_migrated_fortran_library).
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
            target_link_libraries(mumps_common PUBLIC esmumps scotch)
        endif()
        # PT-Scotch C side: mumps_scotch.c's MUMPS_DGRAPHINIT shim calls
        # bare SCOTCH_dgraphInit under -Dptscotch (the define is applied for
        # both languages by the mumps_common ptscotch block above, mirroring
        # -Dscotch), and includes ptscotch.h. Force-include the distributed
        # rename companion (C-only) so that bare call binds to
        # SCOTCH_dgraphInit_mumps. Link the distributed increment before the
        # sequential library.
        if(MUMPS_HAVE_PTSCOTCH)
            target_compile_options(mumps_common PRIVATE
                "SHELL:$<$<COMPILE_LANGUAGE:C>:-include ${_mumps_scotch_inc}/scotch_rename_pt_mumps.h>")
            target_link_libraries(mumps_common PUBLIC ptscotch)
        endif()
        if(MPI_C_FOUND)
            target_link_libraries(mumps_common PRIVATE MPI::MPI_C)
            target_link_libraries(mumps_common INTERFACE
                $<LINK_ONLY:MPI::MPI_C>)
        endif()

        # ── Fold the migrated C bridge INTO the Fortran archive ─────
        # The two per-arithmetic bridge objects are compiled INTO
        # ${LIB_PAIR_PREFIX}mumps, so libeymumps.a carries both the e/y Fortran
        # solvers and their C entry points (${LIB_PREFIX}mumps_c() /
        # ${LIB_PREFIX_COMPLEX}mumps_c()) — no separate libeymumps_c.a.
        # The shared runtime lives inside mumps_common (folded above),
        # which ${LIB_PAIR_PREFIX}mumps already PUBLIC-links — so it is LINKED,
        # not bundled here, and isn't duplicated against the genuine bridge
        # in a combined link. The bridge↔solver reference (qmumps_c →
        # qmumps_f77_) is internal to the one archive; only the
        # solver↔common (runtime) cycle remains, which the
        # ${LIB_PAIR_PREFIX}mumps_full RESCAN below resolves.
        target_sources(${LIB_PAIR_PREFIX}mumps PRIVATE
            $<TARGET_OBJECTS:${LIB_PREFIX}mumps_c_obj>
            $<TARGET_OBJECTS:${LIB_PREFIX_COMPLEX}mumps_c_obj>)
        target_include_directories(${LIB_PAIR_PREFIX}mumps
            PUBLIC
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_target}>
                $<BUILD_INTERFACE:${_mumps_c_bridge_inc_shared}>
                $<BUILD_INTERFACE:${_mumps_c_inc}>
                $<BUILD_INTERFACE:${_mumps_c_src}>
                $<INSTALL_INTERFACE:include/mumps>)
        if(MPI_C_FOUND)
            target_link_libraries(${LIB_PAIR_PREFIX}mumps PRIVATE MPI::MPI_C)
            target_link_libraries(${LIB_PAIR_PREFIX}mumps INTERFACE
                $<LINK_ONLY:MPI::MPI_C>)
        endif()

        # Umbrella INTERFACE target. The qxmumps archive (Fortran solvers
        # + the folded-in C bridge) and mumps_common (Fortran common +
        # the folded-in C runtime) are mutually recursive at link time:
        # qxmumps.a references mumps_size_c_ (in mumps_common's
        # mumps_addr.c) and the common is scanned back against qxmumps.a.
        # A plain left-to-right scan misses one or the other. The
        # ``${LIB_PAIR_PREFIX}mumps_full`` target wraps the archives in a
        # ``LINK_GROUP:RESCAN`` so ld --start-group / --end-group resolves
        # the cycle automatically. Consumers who use
        # ``kind16_libraries::${LIB_PAIR_PREFIX}mumps_full`` don't have to know
        # about the static-archive ordering subtlety. Requires CMake 3.24+
        # on the consumer side. ``mumps_common`` is named explicitly in
        # the group because mumps_size_c_ (Fortran → runtime) lives
        # there, in the folded-in C runtime.
        add_library(${LIB_PAIR_PREFIX}mumps_full INTERFACE)
        # Propagate include directories (mumps_full is a meta-target;
        # by itself it has no headers, so pull them from qxmumps
        # which now carries the per-target bridge include path).
        target_include_directories(${LIB_PAIR_PREFIX}mumps_full INTERFACE
            $<INSTALL_INTERFACE:include/mumps>)
        # ``install(EXPORT)`` rewrites bare target names to
        # ``${PROJECT_NAME}::name`` in INTERFACE_LINK_LIBRARIES, but
        # only for top-level items — it doesn't recurse into the
        # comma-separated arg list of a ``LINK_GROUP`` genex. We
        # therefore embed the namespace literals here so the
        # serialized Targets file points at valid imported targets at
        # consumer time.
        target_link_libraries(${LIB_PAIR_PREFIX}mumps_full INTERFACE
            "$<LINK_GROUP:RESCAN,${PROJECT_NAME}::${LIB_PAIR_PREFIX}mumps,${PROJECT_NAME}::mumps_common>")

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
        # NOTE: genuine s/c/d/z and the extended-precision solvers
        # share ONE common — the migrated ``mumps_common``. Its 36 files
        # are pinned verbatim from upstream (recipe ``copy_files``), so
        # they carry no working precision and serve every arithmetic
        # identically — exactly as upstream compiles ``libmumps_common``
        # once and links it into libsmumps/libdmumps/libcmumps/libzmumps
        # alike. With every common module writing into the single shared
        # Fortran module directory, one producer keeps the .mod set
        # unambiguous for both the genuine and the extended solvers.
        find_package(BLAS QUIET)
        find_package(LAPACK QUIET)
        # ── Parallel stack for the genuine solvers ──────────────────
        # Prefer the in-tree double-precision parallel stack (ScaLAPACK
        # / PBLAS / BLACS / LAPACK / BLAS built above from the SAME MPI
        # this release configured with). Being compiled against the
        # chosen MPI it is ABI-compatible with that MPI by construction,
        # so the genuine build does not hinge on a SYSTEM ScaLAPACK
        # whose BLACS flavor happens to match the active MPI (fragile
        # under, e.g., Intel MPI). Each std archive (plain name, from
        # add_standard_*) is paired with its migrated ``_common``
        # companion (from add_migrated_*); keep whichever candidates are
        # real targets in this configuration.
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
            # as the extended ${LIB_PAIR_PREFIX}mumps does.

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
                # mumps_common / ${LIB_PAIR_PREFIX}mumps).
                target_link_libraries(${_g} PRIVATE MPI::MPI_Fortran)
                # Add_ (C name-mangling) only for the gpu-stub C; don't
                # pass -DAdd_ to the Fortran compile.
                target_compile_definitions(${_g} PRIVATE
                    $<$<COMPILE_LANGUAGE:C>:Add_>)
                # Ordering defines, Fortran-only: the ${ARITH}mumps_gpu.c
                # stubs in these archives never gate on them.
                _add_mumps_ordering_defines(${_g} FORTRAN_ONLY)
                # Pristine upstream compiles within column limits, but
                # relieve line length anyway to match the migrated
                # archives and stay robust to compiler defaults.
                fortran_relax_line_length(${_g})
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
                # PT-Scotch distributed increment (upstream -lptscotch
                # -lscotch order; scotch itself is folded into mumps_common's
                # C runtime). BUILD_INTERFACE so it backs our tests but is
                # not forced into install(EXPORT) — same convention as the
                # migrated ${LIB_PAIR_PREFIX}mumps and the parallel stack below.
                if(MUMPS_HAVE_PTSCOTCH)
                    target_link_libraries(${_s} PUBLIC
                        $<BUILD_INTERFACE:ptscotch>)
                endif()
                if(_have_bundled_parallel)
                    # MPI-matched in-tree stack. The circular refs between
                    # these static archives are resolved by the RESCAN
                    # link group each test builds (tests/mumps); listing
                    # them PUBLIC here mirrors the migrated ${LIB_PAIR_PREFIX}mumps.
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
            foreach(_a d z s c)
                _add_mumps_c_object(${_a}mumps_c_obj ${_a})
            endforeach()

            # Fold the per-family C bridge objects INTO the solver
            # archives: libdzmumps.a carries both the d/z Fortran solvers
            # and their C entry points (dmumps_c()+zmumps_c()), libscmumps.a
            # the s/c pair — no separate libdzmumps_c.a / libscmumps_c.a,
            # matching the folded extended ${LIB_PAIR_PREFIX}mumps. The shared
            # runtime lives inside mumps_common (linked above), not
            # bundled. The bridge↔solver reference is internal to the
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
            # as ${LIB_PAIR_PREFIX}mumps_full; RESCAN over the (combined) solver
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
