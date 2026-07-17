/* multifloats_mpi.cpp -- runtime registration of multifloats MPI handles.
 *
 * Uses the C++ multifloats::float64x2 operators (+, abs via < 0 ? -x : x,
 * comparison) so no manual Knuth two-sum is needed.
 *
 * Compiled as C++; all exported symbols use extern "C" linkage so they
 * are callable from both C and Fortran (through the migrated BLACS
 * wrappers).
 */
#include "multifloats_bridge.h"

extern "C" {

MPI_Datatype MPI_FLOAT64X2    = MPI_DATATYPE_NULL;
MPI_Datatype MPI_COMPLEX64X2 = MPI_DATATYPE_NULL;

MPI_Op MPI_MM_SUM = MPI_OP_NULL;
MPI_Op MPI_WW_SUM = MPI_OP_NULL;
MPI_Op MPI_MM_AMX = MPI_OP_NULL;
MPI_Op MPI_MM_AMN = MPI_OP_NULL;
MPI_Op MPI_WW_AMX = MPI_OP_NULL;
MPI_Op MPI_WW_AMN = MPI_OP_NULL;

/* Fortran-side handles. Populated by multifloats_mpi_init() via
 * MPI_Type_c2f / MPI_Op_c2f after the C-side handles are ready, and
 * surfaced to Fortran via multifloats_mpi_f.f90 using bind(c, name=...).
 * MUMPS calls MPI from Fortran directly with names like
 * ``MPI_FLOAT64X2`` and ``MPI_WW_SUM``, which need INTEGER Fortran
 * handles rather than the C-side MPI_Datatype / MPI_Op opaque types.
 *
 * The two *datatype* handles default to the libmpiseq derived-type
 * sentinel (0x10000000 | total_bytes) rather than 0, so a sequential
 * (libmpiseq) consumer can drive MUMPS with NEITHER MPI_Init NOR
 * multifloats_mpi_init(): the handle already carries the value that
 * libseq's patched MUMPS_COPY dispatches on (float64x2 -> 16 bytes ->
 * 268435472; complex64x2 -> 32 bytes -> 268435488). Without this a
 * skipped init leaves them 0, and any non-in-place reduction STOPs in
 * MUMPS_COPY (DATATYPE=0 matches no branch). In a real-MPI build these
 * defaults are overwritten by multifloats_mpi_init() with the genuine
 * MPI_Type_c2f handles before first use, so the seq default is inert
 * there. The encoding mirrors MPISEQ_DTYPE_TAG in
 * src/mpiseq/mpiseq_c_stubs.c and the MUMPS_COPY cases added by
 * codegen/migrator/libseq_patch.py -- keep the three in sync. */
#define MF_MPISEQ_DTYPE_TAG        0x10000000
#define MF_MPISEQ_SENTINEL(bytes)  (MF_MPISEQ_DTYPE_TAG | (bytes))
MPI_Fint mf_mpi_float64x2_f   = MF_MPISEQ_SENTINEL(16);  /* 268435472 */
MPI_Fint mf_mpi_complex64x2_f = MF_MPISEQ_SENTINEL(32);  /* 268435488 */
MPI_Fint mf_mpi_mm_sum_f      = 0;
MPI_Fint mf_mpi_mm_amx_f      = 0;
MPI_Fint mf_mpi_mm_amn_f      = 0;
MPI_Fint mf_mpi_ww_sum_f      = 0;
MPI_Fint mf_mpi_ww_amx_f      = 0;
MPI_Fint mf_mpi_ww_amn_f      = 0;

} /* extern "C" */

/* ---- User-op callbacks ------------------------------------------ */

/* The static_casts from MPI's void* buffers assume the payload is a
 * contiguous array of float64x2/complex64x2 (plain structs of doubles,
 * 8-byte aligned) — true for every caller, since the buffers originate
 * as Fortran TYPE(real64x2)/TYPE(cmplx64x2) arrays or their C mirrors. */

static void mm_sum_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i) b[i] = b[i] + a[i];
}

static void ww_sum_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i) {
        b[i].re = b[i].re + a[i].re;
        b[i].im = b[i].im + a[i].im;
    }
}

static void mm_amx_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_abs(a[i]) > mf_abs(b[i])) b[i] = a[i];
}

static void mm_amn_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<float64x2 *>(in);
    auto *b = static_cast<float64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_abs(a[i]) < mf_abs(b[i])) b[i] = a[i];
}

static void ww_amx_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_cabs1(a[i]) > mf_cabs1(b[i])) b[i] = a[i];
}

static void ww_amn_fn(void *in, void *inout, int *len, MPI_Datatype *) {
    auto *a = static_cast<complex64x2 *>(in);
    auto *b = static_cast<complex64x2 *>(inout);
    for (int i = 0; i < *len; ++i)
        if (mf_cabs1(a[i]) < mf_cabs1(b[i])) b[i] = a[i];
}

/* ---- One-time registration -------------------------------------- */

/* Under Intel MPI the ops are declared non-commutative even though
 * sum/abs-max/abs-min are commutative: its shm-optimized reduce for
 * commutative user-defined ops hands the callback an 8-byte-aligned
 * pointer into the shm cell payload above the short-message cutoff; a
 * callback whose element type wants 16-byte alignment is compiled with
 * aligned SSE loads and faults on it (#GP -> SIGSEGV, si_addr=0;
 * observed with 2021.18). Non-commutative ops take the sound
 * rank-ordered path through properly aligned buffers. Other MPIs are
 * unaffected (OpenMPI 4.1.6 verified clean) and keep the commutative
 * declaration and its cheaper reduction algorithms. Intel MPI's mpi.h
 * also defines the MPICH_* macros, so key on I_MPI_* only. */
#if defined(I_MPI_VERSION) || defined(I_MPI_NUMVERSION)
#  define EP_MPI_OP_COMMUTE 0
#else
#  define EP_MPI_OP_COMMUTE 1
#endif

extern "C" void multifloats_mpi_init(void) {
    static int initialized = 0;
    if (initialized) return;

    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_FLOAT64X2);
    MPI_Type_commit(&MPI_FLOAT64X2);
    MPI_Type_contiguous(4, MPI_DOUBLE, &MPI_COMPLEX64X2);
    MPI_Type_commit(&MPI_COMPLEX64X2);

    MPI_Op_create(mm_sum_fn, EP_MPI_OP_COMMUTE, &MPI_MM_SUM);
    MPI_Op_create(ww_sum_fn, EP_MPI_OP_COMMUTE, &MPI_WW_SUM);
    MPI_Op_create(mm_amx_fn, EP_MPI_OP_COMMUTE, &MPI_MM_AMX);
    MPI_Op_create(mm_amn_fn, EP_MPI_OP_COMMUTE, &MPI_MM_AMN);
    MPI_Op_create(ww_amx_fn, EP_MPI_OP_COMMUTE, &MPI_WW_AMX);
    MPI_Op_create(ww_amn_fn, EP_MPI_OP_COMMUTE, &MPI_WW_AMN);

    mf_mpi_float64x2_f   = MPI_Type_c2f(MPI_FLOAT64X2);
    mf_mpi_complex64x2_f = MPI_Type_c2f(MPI_COMPLEX64X2);
    mf_mpi_mm_sum_f      = MPI_Op_c2f(MPI_MM_SUM);
    mf_mpi_mm_amx_f      = MPI_Op_c2f(MPI_MM_AMX);
    mf_mpi_mm_amn_f      = MPI_Op_c2f(MPI_MM_AMN);
    mf_mpi_ww_sum_f      = MPI_Op_c2f(MPI_WW_SUM);
    mf_mpi_ww_amx_f      = MPI_Op_c2f(MPI_WW_AMX);
    mf_mpi_ww_amn_f      = MPI_Op_c2f(MPI_WW_AMN);

    initialized = 1;
}
