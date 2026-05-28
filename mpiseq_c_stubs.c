/*
 * C-side MPI stubs for libmpiseq, replacing libseq's bundled mpic.c
 * outright. We compile against the configured MPI's mpi.h (Intel MPI
 * by default, see cmake/CMakeLists.txt) so libmpiseq.a is ABI-
 * compatible with archives compiled against real MPI: same opaque-type
 * sizes, same constant values, same calling conventions. The
 * implementation is single-rank stubs.
 *
 * libseq's mpic.c is dropped from the mpiseq target because Intel's
 * mpi.h defines ``MPI_Comm_f2c`` as a function-like macro; libseq's
 * mpic.c declares it as a function, and the two collide before any
 * compiler tricks can paper over it. This file provides the four
 * routines libseq's mpic.c shipped (MPI_Init, MPI_Comm_rank,
 * MPI_Finalize, MUMPS_CHECKADDREQUAL + aliases, MPI_Wtime) plus every
 * additional MPI_* symbol the BLACS / PBLAS / multifloats_mpi C side
 * calls and which Intel's MPI runtime would normally provide.
 *
 * All stubs assume rank=0, size=1, and either no-op (for collectives,
 * since there's nothing to coordinate with) or print a "should not be
 * called" error and abort (for point-to-point, which a single-rank
 * solver should never invoke). Matches the contract upstream's
 * Fortran-side mpi.f stubs use.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "mpi.h"      /* configured MPI's mpi.h (Intel MPI) */
#include "elapse.h"   /* mumps_elapse() — libseq-supplied wallclock */

/* Intel's mpi.h defines MPI_Comm_f2c as a function-like macro. We need
 * to redefine it as a callable function so Fortran-side mpi.f's
 * ``MPI_Comm_f2c`` reference resolves. Strip the macro before defining
 * the function symbol. */
#ifdef MPI_Comm_f2c
#  undef MPI_Comm_f2c
#endif

static void mpiseq_should_not_be_called(const char *fn)
{
    fprintf(stderr,
            "Error. %s should not be called from a single-rank MUMPS "
            "executable built against libmpiseq.\n", fn);
    abort();
}

/* ── Lifecycle ───────────────────────────────────────────────────────── */

int MPI_Init(int *pargc, char ***pargv)     { (void) pargc; (void) pargv; return MPI_SUCCESS; }
int MPI_Finalize(void)                      { return MPI_SUCCESS; }
int MPI_Comm_rank(MPI_Comm comm, int *rank) { (void) comm; *rank = 0; return MPI_SUCCESS; }
int MPI_Initialized(int *flag)              { *flag = 1; return MPI_SUCCESS; }
int MPI_Comm_size(MPI_Comm comm, int *size) { (void) comm; *size = 1; return MPI_SUCCESS; }
MPI_Comm MPI_Comm_f2c(MPI_Fint comm)        { return (MPI_Comm) comm; }
double MPI_Wtime(void)                      { double v; mumps_elapse(&v); return v; }

/* MUMPS_CHECKADDREQUAL — used by Fortran-side MPI_IS_IN_PLACE check
 * to compare two C addresses, returning 1 if equal. */
void MUMPS_CHECKADDREQUAL(char *a, char *b, int64_t *i)
{
    *i = (a - b == 0) ? 1 : 0;
}
void MUMPS_CHECKADDREQUAL_(char *a, char *b, int64_t *i)
    { MUMPS_CHECKADDREQUAL(a, b, i); }
void mumps_checkaddrequal_(char *a, char *b, int64_t *i)
    { MUMPS_CHECKADDREQUAL(a, b, i); }
void mumps_checkaddrequal__(char *a, char *b, int64_t *i)
    { MUMPS_CHECKADDREQUAL(a, b, i); }
void mumps_checkaddrequal(char *a, char *b, int64_t *i)
    { MUMPS_CHECKADDREQUAL(a, b, i); }

int MPI_Abort(MPI_Comm comm, int errorcode)
{
    (void) comm;
    fprintf(stderr, "MPI_Abort called with errorcode=%d\n", errorcode);
    exit(errorcode);
    return MPI_SUCCESS;
}

/* ── Communicator / group queries ────────────────────────────────────── */

int MPI_Comm_dup(MPI_Comm comm, MPI_Comm *newcomm)         { *newcomm = comm; return MPI_SUCCESS; }
int MPI_Comm_free(MPI_Comm *comm)                          { (void) comm; return MPI_SUCCESS; }
int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm *newcomm)
{
    (void) group; *newcomm = comm; return MPI_SUCCESS;
}
int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm *newcomm)
{
    (void) color; (void) key; *newcomm = comm; return MPI_SUCCESS;
}
int MPI_Comm_group(MPI_Comm comm, MPI_Group *group)        { (void) comm; *group = 0; return MPI_SUCCESS; }
int MPI_Comm_get_attr(MPI_Comm c, int kv, void *v, int *f) { (void) c; (void) kv; (void) v; *f = 0; return MPI_SUCCESS; }
int MPI_Group_incl(MPI_Group g, int n, const int r[], MPI_Group *ng)
{
    (void) g; (void) n; (void) r; *ng = 0; return MPI_SUCCESS;
}
int MPI_Group_free(MPI_Group *g)                           { *g = 0; return MPI_SUCCESS; }

/* ── Collectives — single rank, so the "operation" is a no-op ────────── */

int MPI_Barrier(MPI_Comm comm)                             { (void) comm; return MPI_SUCCESS; }

int MPI_Bcast(void *buf, int count, MPI_Datatype dt, int root, MPI_Comm comm)
{
    (void) buf; (void) count; (void) dt; (void) root; (void) comm; return MPI_SUCCESS;
}

/* Single-rank reduce semantically equals ``recvbuf := sendbuf``. With
 * the configured MPI's mpi.h on the include path the constants and
 * opaque types match what callers see, but we still don't carry a
 * datatype-size table here — element size for the memcpy isn't known
 * without a dispatch on every Intel MPI MPI_Datatype value. Abort on
 * any non-Fortran call; the Fortran-side mpi.f handles the in-place
 * + dispatch case where MPI_IN_PLACE / MPI_DOUBLE_PRECISION / etc.
 * COMMON-block constants are visible. */
int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,
               MPI_Datatype dt, MPI_Op op, int root, MPI_Comm comm)
{
    (void) sendbuf; (void) recvbuf; (void) count;
    (void) dt; (void) op; (void) root; (void) comm;
    mpiseq_should_not_be_called("MPI_Reduce");
    return MPI_SUCCESS;
}

int MPI_Allreduce(const void *sendbuf, void *recvbuf, int count,
                  MPI_Datatype dt, MPI_Op op, MPI_Comm comm)
{
    (void) sendbuf; (void) recvbuf; (void) count;
    (void) dt; (void) op; (void) comm;
    mpiseq_should_not_be_called("MPI_Allreduce");
    return MPI_SUCCESS;
}

/* ── Point-to-point — illegal for a single rank ──────────────────────── */

int MPI_Send(const void *b, int c, MPI_Datatype d, int dst, int tag, MPI_Comm cm)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    mpiseq_should_not_be_called("MPI_Send");
    return MPI_SUCCESS;
}
int MPI_Rsend(const void *b, int c, MPI_Datatype d, int dst, int tag, MPI_Comm cm)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    mpiseq_should_not_be_called("MPI_Rsend");
    return MPI_SUCCESS;
}
int MPI_Recv(void *b, int c, MPI_Datatype d, int src, int tag,
             MPI_Comm cm, MPI_Status *st)
{
    (void) b; (void) c; (void) d; (void) src; (void) tag; (void) cm; (void) st;
    mpiseq_should_not_be_called("MPI_Recv");
    return MPI_SUCCESS;
}
int MPI_Isend(const void *b, int c, MPI_Datatype d, int dst, int tag,
              MPI_Comm cm, MPI_Request *rq)
{
    (void) b; (void) c; (void) d; (void) dst; (void) tag; (void) cm;
    if (rq) *rq = 0;
    mpiseq_should_not_be_called("MPI_Isend");
    return MPI_SUCCESS;
}
int MPI_Irecv(void *b, int c, MPI_Datatype d, int src, int tag,
              MPI_Comm cm, MPI_Request *rq)
{
    (void) b; (void) c; (void) d; (void) src; (void) tag; (void) cm;
    if (rq) *rq = 0;
    mpiseq_should_not_be_called("MPI_Irecv");
    return MPI_SUCCESS;
}
int MPI_Sendrecv(const void *sb, int sc, MPI_Datatype sd, int dst, int stag,
                 void *rb, int rc, MPI_Datatype rd, int src, int rtag,
                 MPI_Comm cm, MPI_Status *st)
{
    (void) sb; (void) sc; (void) sd; (void) dst; (void) stag;
    (void) rb; (void) rc; (void) rd; (void) src; (void) rtag;
    (void) cm; (void) st;
    mpiseq_should_not_be_called("MPI_Sendrecv");
    return MPI_SUCCESS;
}
int MPI_Waitall(int count, MPI_Request reqs[], MPI_Status stats[])
{
    (void) count; (void) reqs; (void) stats; return MPI_SUCCESS;
}
int MPI_Testall(int count, MPI_Request reqs[], int *flag, MPI_Status stats[])
{
    (void) count; (void) reqs; (void) stats; *flag = 1; return MPI_SUCCESS;
}

/* ── Datatype handles ────────────────────────────────────────────────── */

/* Return element-size of an Intel-MPI base datatype constant, or 0 if
 * not a recognised base type. The Intel constants encode size in bits
 * 16-23 (e.g. MPI_DOUBLE = 0x4c00080b → 0x08 = 8 bytes), but we keep
 * the dispatch explicit for clarity and to cover the named-Fortran-
 * variant constants whose IDs use a different layout. */
static int mpiseq_base_type_bytes(MPI_Datatype t)
{
    if (t == MPI_BYTE)             return 1;
    if (t == MPI_CHAR)             return 1;
    if (t == MPI_SHORT)            return 2;
    if (t == MPI_INT)              return 4;
    if (t == MPI_INTEGER)          return 4;
    if (t == MPI_LONG)             return sizeof(long);
    if (t == MPI_LONG_LONG_INT)    return 8;
    if (t == MPI_INTEGER8)         return 8;
    if (t == MPI_FLOAT)            return 4;
    if (t == MPI_REAL)             return 4;
    if (t == MPI_DOUBLE)           return 8;
    if (t == MPI_DOUBLE_PRECISION) return 8;
    if (t == MPI_REAL8)            return 8;
    if (t == MPI_LONG_DOUBLE)      return 10;   /* x87 80-bit */
    if (t == MPI_REAL16)           return 16;
    if (t == MPI_COMPLEX)          return 8;
    if (t == MPI_DOUBLE_COMPLEX)   return 16;
    if (t == MPI_COMPLEX32)        return 32;
    if (t == MPI_C_LONG_DOUBLE_COMPLEX) return 20;
    /* Already-derived multifloats sentinels: pass byte-size through. */
    if (((uintptr_t) t & 0xFF000000U) == 0x10000000U)
        return (int)((uintptr_t) t & 0x00FFFFFFU);
    return 0;
}

/* libmpiseq derived-datatype sentinels: 0x10000000 | total_size_in_bytes.
 * MUMPS_COPY in libseq's mpi.f recognises this tag and dispatches on
 * the encoded size, so multifloats' MPI_FLOAT64X2 / MPI_COMPLEX64X2
 * survive a single-rank MPI_ALLREDUCE without the C-side user-op
 * callbacks ever firing. The 0x10000000 tag is below MPI_OP_NULL
 * (0x18000000) and well clear of Intel's 0x4c00**** datatype range. */
#define MPISEQ_DTYPE_TAG  0x10000000U

int MPI_Type_contiguous(int count, MPI_Datatype old, MPI_Datatype *new_)
{
    int b = mpiseq_base_type_bytes(old);
    *new_ = (MPI_Datatype)(uintptr_t)(MPISEQ_DTYPE_TAG | (unsigned)(count * b));
    return MPI_SUCCESS;
}
int MPI_Type_vector(int c, int bl, int st, MPI_Datatype old, MPI_Datatype *new_)
{
    (void) st;
    int b = mpiseq_base_type_bytes(old);
    *new_ = (MPI_Datatype)(uintptr_t)(MPISEQ_DTYPE_TAG | (unsigned)(c * bl * b));
    return MPI_SUCCESS;
}
int MPI_Type_indexed(int c, const int bls[], const int disps[],
                     MPI_Datatype old, MPI_Datatype *new_)
{
    (void) disps;
    int b = mpiseq_base_type_bytes(old);
    int total = 0;
    for (int i = 0; i < c; ++i) total += bls[i];
    *new_ = (MPI_Datatype)(uintptr_t)(MPISEQ_DTYPE_TAG | (unsigned)(total * b));
    return MPI_SUCCESS;
}
int MPI_Type_create_struct(int c, const int bls[], const MPI_Aint disps[],
                           const MPI_Datatype types[], MPI_Datatype *new_)
{
    (void) disps;
    int total = 0;
    for (int i = 0; i < c; ++i) total += bls[i] * mpiseq_base_type_bytes(types[i]);
    *new_ = (MPI_Datatype)(uintptr_t)(MPISEQ_DTYPE_TAG | (unsigned)total);
    return MPI_SUCCESS;
}
int MPI_Type_commit(MPI_Datatype *t)                       { (void) t; return MPI_SUCCESS; }
int MPI_Type_free(MPI_Datatype *t)                         { (void) t; return MPI_SUCCESS; }
int MPI_Type_match_size(int typeclass, int size, MPI_Datatype *t)
{
    (void) typeclass; (void) size; *t = MPI_DATATYPE_NULL; return MPI_SUCCESS;
}

int MPI_Pack_size(int incount, MPI_Datatype dt, MPI_Comm cm, int *size)
{
    (void) dt; (void) cm; *size = incount; return MPI_SUCCESS;
}
int MPI_Pack(const void *inbuf, int incount, MPI_Datatype dt,
             void *outbuf, int outsize, int *position, MPI_Comm cm)
{
    (void) dt; (void) outsize; (void) cm;
    if (inbuf && outbuf && incount > 0)
        memcpy((char *) outbuf + *position, inbuf, (size_t) incount);
    if (position) *position += incount;
    return MPI_SUCCESS;
}
int MPI_Unpack(const void *inbuf, int insize, int *position,
               void *outbuf, int outcount, MPI_Datatype dt, MPI_Comm cm)
{
    (void) insize; (void) dt; (void) cm;
    if (inbuf && outbuf && outcount > 0)
        memcpy(outbuf, (const char *) inbuf + *position, (size_t) outcount);
    if (position) *position += outcount;
    return MPI_SUCCESS;
}

/* libmpiseq op sentinels: distinct non-null handles above MPI_OP_NULL.
 * We never actually invoke the user callback (single-rank ALLREDUCE
 * collapses to memcpy in MUMPS_COPY); callers just need a non-null
 * handle to round-trip through MPI_Allreduce. */
int MPI_Op_create(MPI_User_function *fn, int commute, MPI_Op *op)
{
    (void) fn; (void) commute;
    static unsigned next = 1;
    *op = (MPI_Op)(uintptr_t)(0x18000000U + next++);
    return MPI_SUCCESS;
}
int MPI_Op_free(MPI_Op *op)                                { *op = MPI_OP_NULL; return MPI_SUCCESS; }

int MPI_Error_class(int errcode, int *errclass)            { *errclass = errcode; return MPI_SUCCESS; }
