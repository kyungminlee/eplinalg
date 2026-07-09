# libmpiseq — sequential MPI stub from MUMPS 5.8.2 upstream (libseq),
# extended with first-party stubs for the extended-precision stack.
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
# The upstream libseq sources are staged into ``_mpiseq_src/`` at the
# staging root by ``migrator stage`` (see ``_std_dirs`` in
# ``src/migrator/staging.py``); the parent gates this subdirectory on
# that directory's presence plus the Intel-MPI ABI check
# (MPISEQ_IS_INTEL_MPI — the stubs' constants are Intel-MPI-specific).

set(_mpiseq_upstream ${CMAKE_SOURCE_DIR}/_mpiseq_src)

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
    file(REMOVE ${_mpiseq_upstream}/mpif.h)
endif()
if(MPI_C_FOUND)
    file(REMOVE ${_mpiseq_upstream}/mpi.h)
endif()
add_library(mpiseq STATIC
    ${_mpiseq_upstream}/mpi.f
    # libseq's bundled mpic.c is replaced by mpiseq_c_stubs.c (added
    # below). Intel MPI defines MPI_Comm_f2c as a function-like macro;
    # libseq's mpic.c declares it as a function and the two collide
    # before any compile-flag tricks. Replacing the whole file is
    # cleaner.
    ${_mpiseq_upstream}/elapse.c
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
# PT-Scotch dgraph stand-ins. When PT-Scotch is staged the mumps
# archives are compiled -Dptscotch and reference the distributed
# SCOTCHFDGRAPH*_MUMPS / SCOTCH_dgraphInit_mumps entry points, which a
# libmpiseq link can't take from the real ptscotch archive (it needs a
# real MPI underneath). Single-rank execution never reaches them, so
# aborting stubs close the link. Always folded in like the C stubs:
# without PT-Scotch nothing references them and the archive members
# stay unpulled (sequential scotch defines no dgraph symbols, so there
# is no overlap either way).
target_sources(mpiseq PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mpiseq_ptscotch_stubs.f90)
# ep_-privatized twins of mpi.f's BLACS / DESCINIT abort stubs. The
# migrated extended-precision MUMPS archives bind the privatized
# ep_ BLACS engine names (recipes/mumps.yaml extra_renames), so a
# libmpiseq link must supply the same never-executed abort stubs
# under the ep_ spelling. Always folded in like the C stubs:
# baseline kind4/kind8 consumers never reference ep_* and the
# archive member stays unpulled.
target_sources(mpiseq PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/mpiseq_ep_stubs.f)
target_include_directories(mpiseq PUBLIC
    $<BUILD_INTERFACE:${_mpiseq_upstream}>
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

# The libmpiseq (seq) release ships libmpiseq so consumers of the
# ``-seq``-tagged MUMPS archives have a sequential MPI stub to link.
# No flavor tag: there is only ever one mpiseq (no real-MPI twin to
# collide with), so it installs as libmpiseq-${FORTRAN_COMPILER_TAG}.a.
# Guarded by the option, so the default in-tree-test-only build of
# mpiseq is unchanged. (Epilogue stays explicit: this subdirectory
# executes before the parent defines its install helpers.)
if(MUMPS_LIBSEQ_RELEASE)
    fortran_install_library(mpiseq
        EXPORT mpiseqTargets NAMESPACE eplinalg::)
    get_target_property(_mpiseq_moddir mpiseq Fortran_MODULE_DIRECTORY)
    if(_mpiseq_moddir)
        fortran_install_modules(mpiseq)
    endif()
endif()
