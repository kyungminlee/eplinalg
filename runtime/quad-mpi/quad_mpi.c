/* quad_mpi.c -- runtime registration of the custom MPI reduction ops
 * used by the migrated kind16 (REAL(KIND=16) / __float128) MUMPS.
 *
 * Unlike the multifloats bridge, the *datatypes* here are the standard,
 * predefined MPI handles MPI_REAL16 / MPI_COMPLEX32 -- mpif.h declares
 * them and every MPI transports them correctly for point-to-point and
 * broadcast. What is NOT safe is a *reduction* on them: Intel MPI 2021's
 * builtin MPI_SUM / MPI_MAX / MPI_MIN have no 16-byte-real kernel and
 * segfault the moment MUMPS issues MPI_REDUCE(..., MPI_REAL16, MPI_SUM,
 * ...) at np >= 2. So we keep the datatypes but register user-defined
 * combine ops via MPI_Op_create and route every migrated reduction
 * through them -- mirroring how multifloats uses MPI_DD_SUM / MPI_DD_AMX.
 *
 * Compiled as plain C (kind16 BLACS/PBLAS also compile as C, using the
 * __float128 native operators). The exported symbols are ordinary C
 * externs; the Fortran side (quad_mpi_f.f90) binds to the *_f handles
 * via bind(c, name=...).
 */
#include <mpi.h>
#include <quadmath.h>
#include <string.h>

typedef __float128 qreal;
typedef struct { __float128 re, im; } xcomplex;  /* == MPI_COMPLEX32 layout */

/* ---- C-side op handles (opaque MPI_Op) --------------------------- */

MPI_Op MPI_QQ_SUM = MPI_OP_NULL;
MPI_Op MPI_QQ_AMX = MPI_OP_NULL;
MPI_Op MPI_QQ_AMN = MPI_OP_NULL;
MPI_Op MPI_XX_SUM = MPI_OP_NULL;
MPI_Op MPI_XX_AMX = MPI_OP_NULL;
MPI_Op MPI_XX_AMN = MPI_OP_NULL;

/* ---- Fortran-side handles ---------------------------------------- *
 * Populated by quad_mpi_init() via MPI_Op_c2f after the C-side ops
 * are created, and surfaced to Fortran via quad_mpi_f.f90 using
 * bind(c, name=...). MUMPS calls MPI from Fortran directly with names
 * like MPI_QQ_SUM, which need INTEGER Fortran handles rather than the
 * C-side opaque MPI_Op. Values are 0 until quad_mpi_init runs. */
MPI_Fint q_mpi_qq_sum_f = 0;
MPI_Fint q_mpi_qq_amx_f = 0;
MPI_Fint q_mpi_qq_amn_f = 0;
MPI_Fint q_mpi_xx_sum_f = 0;
MPI_Fint q_mpi_xx_amx_f = 0;
MPI_Fint q_mpi_xx_amn_f = 0;

/* ---- User-op callbacks ------------------------------------------- *
 * MUMPS reduces non-negative norms / scaling factors, so abs-max and
 * abs-min coincide with true max/min there; taking the magnitude keeps
 * the op well defined for any input and matches the multifloats sibling
 * (MPI_DD_AMX / MPI_DD_AMN).
 *
 * ALIGNMENT: __float128 arithmetic uses SSE aligned moves (movaps),
 * which fault on operands that are not 16-byte aligned. MPI hands the
 * user-op scratch buffers whose alignment is NOT guaranteed to match the
 * type -- Intel MPI's intranode "tiny" reduce path in particular passes
 * an 8-byte-aligned staging buffer for MPI_REAL16. Dereferencing it
 * directly as ``qreal *`` SIGSEGVs at np >= 2 (np == 1 / MPISEQ never
 * invoke the op, so they are unaffected). We therefore memcpy every
 * element through 16-byte-aligned stack locals -- the compiler emits an
 * unaligned load/store for the memcpy but keeps the arithmetic on the
 * aligned automatics. */

static __float128 xcabs1(xcomplex z) { return fabsq(z.re) + fabsq(z.im); }

static void qq_sum_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(qreal), pb += sizeof(qreal)) {
        qreal a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        b = b + a;
        memcpy(pb, &b, sizeof b);
    }
}

static void qq_amx_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(qreal), pb += sizeof(qreal)) {
        qreal a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        if (fabsq(a) > fabsq(b)) memcpy(pb, &a, sizeof a);
    }
}

static void qq_amn_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(qreal), pb += sizeof(qreal)) {
        qreal a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        if (fabsq(a) < fabsq(b)) memcpy(pb, &a, sizeof a);
    }
}

static void xx_sum_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(xcomplex), pb += sizeof(xcomplex)) {
        xcomplex a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        b.re = b.re + a.re;
        b.im = b.im + a.im;
        memcpy(pb, &b, sizeof b);
    }
}

static void xx_amx_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(xcomplex), pb += sizeof(xcomplex)) {
        xcomplex a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        if (xcabs1(a) > xcabs1(b)) memcpy(pb, &a, sizeof a);
    }
}

static void xx_amn_fn(void *in, void *inout, int *len, MPI_Datatype *dt) {
    char *pa = (char *)in, *pb = (char *)inout;
    (void)dt;
    for (int i = 0; i < *len; ++i, pa += sizeof(xcomplex), pb += sizeof(xcomplex)) {
        xcomplex a, b;
        memcpy(&a, pa, sizeof a);
        memcpy(&b, pb, sizeof b);
        if (xcabs1(a) < xcabs1(b)) memcpy(pb, &a, sizeof a);
    }
}

/* ---- One-time registration --------------------------------------- */

/* Under Intel MPI the ops are declared non-commutative even though
 * sum/abs-max/abs-min are commutative: its shm-optimized reduce for
 * commutative user-defined ops hands the callback an 8-byte-aligned
 * pointer into the shm cell payload above the short-message cutoff;
 * these ops load __float128 with aligned SSE instructions and fault on
 * it (#GP -> SIGSEGV, si_addr=0; observed with 2021.18). Non-commutative
 * ops take the sound rank-ordered path through properly aligned buffers.
 * Other MPIs are unaffected (OpenMPI 4.1.6 verified clean) and keep the
 * commutative declaration and its cheaper reduction algorithms. Intel
 * MPI's mpi.h also defines the MPICH_* macros, so key on I_MPI_* only. */
#if defined(I_MPI_VERSION) || defined(I_MPI_NUMVERSION)
#  define EP_MPI_OP_COMMUTE 0
#else
#  define EP_MPI_OP_COMMUTE 1
#endif

void quad_mpi_init(void) {
    static int initialized = 0;
    if (initialized) return;

    MPI_Op_create(qq_sum_fn, EP_MPI_OP_COMMUTE, &MPI_QQ_SUM);
    MPI_Op_create(qq_amx_fn, EP_MPI_OP_COMMUTE, &MPI_QQ_AMX);
    MPI_Op_create(qq_amn_fn, EP_MPI_OP_COMMUTE, &MPI_QQ_AMN);
    MPI_Op_create(xx_sum_fn, EP_MPI_OP_COMMUTE, &MPI_XX_SUM);
    MPI_Op_create(xx_amx_fn, EP_MPI_OP_COMMUTE, &MPI_XX_AMX);
    MPI_Op_create(xx_amn_fn, EP_MPI_OP_COMMUTE, &MPI_XX_AMN);

    q_mpi_qq_sum_f = MPI_Op_c2f(MPI_QQ_SUM);
    q_mpi_qq_amx_f = MPI_Op_c2f(MPI_QQ_AMX);
    q_mpi_qq_amn_f = MPI_Op_c2f(MPI_QQ_AMN);
    q_mpi_xx_sum_f = MPI_Op_c2f(MPI_XX_SUM);
    q_mpi_xx_amx_f = MPI_Op_c2f(MPI_XX_AMX);
    q_mpi_xx_amn_f = MPI_Op_c2f(MPI_XX_AMN);

    initialized = 1;
}
