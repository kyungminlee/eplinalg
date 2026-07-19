/* mpiseq_dtype_tag.h -- libmpiseq derived-datatype sentinel encoding.
 *
 * libmpiseq derived-datatype sentinels: 0x10000000 | total_size_in_bytes.
 * MUMPS_COPY in libseq's mpi.f recognises this tag and dispatches on
 * the encoded size, so multifloats' MPI_FLOAT64X2 / MPI_COMPLEX64X2
 * survive a single-rank MPI_ALLREDUCE without the C-side user-op
 * callbacks ever firing. The 0x10000000 tag is below MPI_OP_NULL
 * (0x18000000) and well clear of Intel's 0x4c00**** datatype range.
 *
 * Shared by mpiseq_c_stubs.c (encodes the sentinels in the MPI_Type_*
 * stubs) and multifloats_mpi.cpp (seq-safe defaults for the Fortran
 * datatype handles). The MUMPS_COPY cases emitted by
 * codegen/migrator/libseq_patch.py bake the decimal expansions
 * (float64x2 -> 16 bytes -> 268435472; complex64x2 -> 32 bytes ->
 * 268435488) into generated Fortran -- keep that file in sync with
 * this header. */
#ifndef MPISEQ_DTYPE_TAG_H
#define MPISEQ_DTYPE_TAG_H

#define MPISEQ_DTYPE_TAG             0x10000000U
#define MPISEQ_DTYPE_SENTINEL(bytes) (MPISEQ_DTYPE_TAG | (unsigned)(bytes))

#endif /* MPISEQ_DTYPE_TAG_H */
