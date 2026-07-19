/* ep_mpi_op_commute.h -- commute flag for the extended-precision
 * user-defined MPI reduction ops. Both runtime bridges (quad-mpi's
 * MPI_QQ_* / MPI_XX_* and multifloats-mpi's MPI_MM_* / MPI_WW_*) pass
 * EP_MPI_OP_COMMUTE to every MPI_Op_create.
 *
 * Under Intel MPI the ops are declared non-commutative even though
 * sum/abs-max/abs-min are commutative: its shm-optimized reduce for
 * commutative user-defined ops hands the callback an 8-byte-aligned
 * pointer into the shm cell payload above the short-message cutoff; a
 * callback whose element type wants 16-byte alignment is compiled with
 * aligned SSE loads and faults on it (#GP -> SIGSEGV, si_addr=0;
 * observed with 2021.18). Non-commutative ops take the sound
 * rank-ordered path through properly aligned buffers. Other MPIs are
 * unaffected (OpenMPI 4.1.6 verified clean) and keep the commutative
 * declaration and its cheaper reduction algorithms. Intel MPI's mpi.h
 * also defines the MPICH_* macros, so key on I_MPI_* only.
 *
 * Include AFTER <mpi.h>: the I_MPI_* detection macros come from Intel's
 * mpi.h. Lives under mpiseq/ because that is the one runtime directory
 * the migrator stages for every target; the two bridge directories are
 * staged per-arithmetic, so neither can host a header the other
 * includes. */
#ifndef EP_MPI_OP_COMMUTE_H
#define EP_MPI_OP_COMMUTE_H

#ifndef MPI_VERSION
#  error "include <mpi.h> before ep_mpi_op_commute.h"
#endif

#if defined(I_MPI_VERSION) || defined(I_MPI_NUMVERSION)
#  define EP_MPI_OP_COMMUTE 0
#else
#  define EP_MPI_OP_COMMUTE 1
#endif

#endif /* EP_MPI_OP_COMMUTE_H */
