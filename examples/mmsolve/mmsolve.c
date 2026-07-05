/*
 * mmsolve — a small MPI sparse direct solver built on MUMPS, linking ALL TEN
 * arithmetics into one binary with runtime `-t` selection.
 *
 *   mmsolve -t <s|c|d|z|e|y|q|x|m|w> [-v] <matrix.mtx> <rhs.mtx> <solution.mtx>
 *
 * Reads a sparse matrix A (Matrix Market coordinate) and a right-hand side b
 * (Matrix Market array, n x 1), solves A x = b with MUMPS, and writes x as a
 * Matrix Market array file. The `-t` type prefix selects the arithmetic:
 *
 *     s = single real          c = single complex
 *     d = double real          z = double complex
 *     e = long-double real     y = long-double complex     (kind10)
 *     q = __float128 real      x = __float128 complex       (kind16)
 *     m = double-double real   w = double-double complex    (multifloats)
 *
 * s/c/d/z are the four "genuine" arithmetics MUMPS ships and the four Intel MKL
 * provides a ScaLAPACK/BLACS backend for. e/y/q/x/m/w are eplinalg's extended
 * precisions. This is the whole point of the Task 29 archive split: the migrator
 * fully prefix-renames every extended stack (emumps/ymumps/qmumps/xmumps/mmumps/
 * wmumps_c, ey/qx/mw-prefixed ScaLAPACK/LAPACK/BLAS/BLACS) and all share ONE
 * arith-agnostic mumps_common, so the ten stacks carry pairwise-disjoint symbols
 * and coexist in a single link with MKL — no -Wl,--allow-multiple-definition.
 * Genuine s/c/d/z run on MKL; the six extended run on their in-tree reference
 * stacks; the type-agnostic plumbing is captured by MKL (linked first).
 *
 * MPI: centralized input. Rank 0 reads the files and owns the assembled
 * matrix + RHS + solution; every rank calls MUMPS collectively. Run with e.g.
 *   mpirun -n 4 ./mmsolve -t d A.mtx b.mtx x.mtx
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

/* Genuine MUMPS C headers (installed under include/mumps by eplinalg). Each
 * declares its own <?>MUMPS_STRUC_C struct and <?>mumps_c() entry point; all
 * four share the one mumps_common runtime, so they coexist in one program. */
#include "smumps_c.h"
#include "cmumps_c.h"
#include "dmumps_c.h"
#include "zmumps_c.h"

/* Extended kind10 arithmetics (long double). Unlike the genuine four, these
 * come from eplinalg's migrated MUMPS stack — but the migrator fully prefix-
 * renames them (emumps_c/ymumps_c, ey-prefixed ScaLAPACK/LAPACK/BLAS) and they
 * share the same arith-agnostic mumps_common, so they carry ZERO symbols that
 * collide with the genuine s/c/d/z and link into this same binary. e = real
 * long double, y = complex long double. Their typedefs (mumps_long_double,
 * mumps_long_double_complex) are guarded, so co-including both is safe. */
#include "emumps_c.h"
#include "ymumps_c.h"

/* Extended kind16 (__float128) and multifloats (double-double) arithmetics.
 * Same story: fully prefix-renamed stacks sharing the one mumps_common. Each
 * pair's shared value typedefs are include-guarded (MUMPS_FLOAT128_TYPES,
 * MUMPS_FLOAT64X2_TYPES), so all ten headers co-include without redefinition.
 * q = __float128 real, x = __float128 complex; m = double-double real,
 * w = double-double complex (real64x2 is a {double limbs[2]} two-limb struct). */
#include "qmumps_c.h"
#include "xmumps_c.h"
#include "mmumps_c.h"
#include "wmumps_c.h"

#include "mmio_min.h"

/* Extended-precision custom MPI reductions. e/y (long double) reduce through
 * native MPI_LONG_DOUBLE, but q/x (__float128) and m/w (real64x2) use custom
 * MPI_Op/datatypes that must be registered after MPI_Init and before the first
 * MUMPS call. The extended BLACS's blacs_pinfo_ normally does this on first
 * use, but in this MKL-first link MKL's blacs_pinfo_ wins, so we register them
 * ourselves. Both are idempotent and safe to call unconditionally. */
extern void quad_mpi_init(void);         /* kind16   q/x custom quad reduce ops */
extern void multifloats_mpi_init(void);  /* multifloats m/w custom reduce ops   */

/* MUMPS job constants and the "no output" ICNTL settings. */
enum { JOB_INIT = -1, JOB_END = -2, JOB_ANALYZE_FACTOR_SOLVE = 6 };

static void quiet_icntl(int *icntl) {
    /* ICNTL(1..3) output streams off, ICNTL(4)=0 verbosity. C is 0-based. */
    icntl[0] = -1; icntl[1] = -1; icntl[2] = -1; icntl[3] = 0;
}

/* ── real solve (s, d) ────────────────────────────────────────────────
 * On the host, assemble typed COO arrays from the double-valued MM, run the
 * MUMPS init→solve→finalize cycle, and copy the centralized solution (which
 * MUMPS returns in id.rhs) back into xr. `sym` (0 or 2) is passed uniformly
 * to every rank, since it must agree across the communicator at analysis. */
#define GEN_REAL(NAME, STRUC, CFUN, REAL)                                        \
static int NAME(const MM *A, const MM *b, int is_host, int sym, int verbose,     \
                double *xr, double *xi) {                                        \
    (void)xi;                                                                    \
    STRUC id; memset(&id, 0, sizeof id);                                         \
    MUMPS_INT *irn = NULL, *jcn = NULL; REAL *a = NULL, *rhs = NULL;             \
    int n = 0; long nz = 0; int infog1 = 0;                                      \
    if (is_host) {                                                               \
        n = A->n; nz = A->nnz;                                                   \
        irn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        jcn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        a   = malloc(sizeof(REAL) * (size_t)nz);                                 \
        rhs = malloc(sizeof(REAL) * (size_t)n);                                  \
        for (long k = 0; k < nz; k++) { irn[k] = A->irn[k]; jcn[k] = A->jcn[k];  \
                                        a[k] = (REAL)A->val[k]; }                \
        for (int i = 0; i < n; i++) rhs[i] = (REAL)b->val[i];                    \
    }                                                                            \
    id.par = 1; id.sym = sym;                                                    \
    id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(MPI_COMM_WORLD);                   \
    id.job = JOB_INIT; CFUN(&id);                                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; goto cleanup; }                 \
    if (!verbose) quiet_icntl(id.icntl);                                         \
    if (is_host) { id.n = n; id.nnz = (MUMPS_INT8)nz;                            \
                   id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs; }         \
    id.job = JOB_ANALYZE_FACTOR_SOLVE; CFUN(&id);                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; id.job = JOB_END; CFUN(&id); goto cleanup; } \
    if (is_host) for (int i = 0; i < n; i++) xr[i] = (double)rhs[i];             \
    id.job = JOB_END; CFUN(&id);                                                 \
cleanup:                                                                         \
    free(irn); free(jcn); free(a); free(rhs);                                    \
    return infog1;                                                               \
}

/* ── complex solve (c, z) ─────────────────────────────────────────────
 * CPLX is the interleaved {r,i} struct MUMPS uses (mumps_complex for c,
 * mumps_double_complex for z). A real MM matrix is promoted with zero
 * imaginary part. */
#define GEN_CPLX(NAME, STRUC, CFUN, CPLX, REAL)                                  \
static int NAME(const MM *A, const MM *b, int is_host, int sym, int verbose,     \
                double *xr, double *xi) {                                        \
    STRUC id; memset(&id, 0, sizeof id);                                         \
    MUMPS_INT *irn = NULL, *jcn = NULL; CPLX *a = NULL, *rhs = NULL;             \
    int n = 0; long nz = 0; int infog1 = 0;                                      \
    if (is_host) {                                                               \
        n = A->n; nz = A->nnz;                                                   \
        irn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        jcn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        a   = malloc(sizeof(CPLX) * (size_t)nz);                                 \
        rhs = malloc(sizeof(CPLX) * (size_t)n);                                  \
        for (long k = 0; k < nz; k++) { irn[k] = A->irn[k]; jcn[k] = A->jcn[k];  \
            a[k].r = (REAL)A->val[k]; a[k].i = (REAL)(A->ival ? A->ival[k] : 0.0); } \
        for (int i = 0; i < n; i++) {                                            \
            rhs[i].r = (REAL)b->val[i]; rhs[i].i = (REAL)(b->ival ? b->ival[i] : 0.0); } \
    }                                                                            \
    id.par = 1; id.sym = sym;                                                    \
    id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(MPI_COMM_WORLD);                   \
    id.job = JOB_INIT; CFUN(&id);                                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; goto cleanup; }                 \
    if (!verbose) quiet_icntl(id.icntl);                                         \
    if (is_host) { id.n = n; id.nnz = (MUMPS_INT8)nz;                            \
                   id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs; }         \
    id.job = JOB_ANALYZE_FACTOR_SOLVE; CFUN(&id);                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; id.job = JOB_END; CFUN(&id); goto cleanup; } \
    if (is_host) for (int i = 0; i < n; i++) { xr[i] = (double)rhs[i].r; xi[i] = (double)rhs[i].i; } \
    id.job = JOB_END; CFUN(&id);                                                 \
cleanup:                                                                         \
    free(irn); free(jcn); free(a); free(rhs);                                    \
    return infog1;                                                               \
}

/* ── multifloats double-double solves (m, w) ──────────────────────────
 * real64x2 is a two-limb {hi,lo} struct, so the scalar cast in GEN_REAL/
 * GEN_CPLX can't build it. Promote each double d to the exact double-double
 * {d, 0} (hi=d, lo=0) on the way in, and read the high limb back out. */
#define GEN_REAL_DD(NAME, STRUC, CFUN, REAL)                                     \
static int NAME(const MM *A, const MM *b, int is_host, int sym, int verbose,     \
                double *xr, double *xi) {                                        \
    (void)xi;                                                                    \
    STRUC id; memset(&id, 0, sizeof id);                                         \
    MUMPS_INT *irn = NULL, *jcn = NULL; REAL *a = NULL, *rhs = NULL;             \
    int n = 0; long nz = 0; int infog1 = 0;                                      \
    if (is_host) {                                                               \
        n = A->n; nz = A->nnz;                                                   \
        irn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        jcn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        a   = malloc(sizeof(REAL) * (size_t)nz);                                 \
        rhs = malloc(sizeof(REAL) * (size_t)n);                                  \
        for (long k = 0; k < nz; k++) { irn[k] = A->irn[k]; jcn[k] = A->jcn[k];  \
            a[k].limbs[0] = A->val[k]; a[k].limbs[1] = 0.0; }                    \
        for (int i = 0; i < n; i++) { rhs[i].limbs[0] = b->val[i]; rhs[i].limbs[1] = 0.0; } \
    }                                                                            \
    id.par = 1; id.sym = sym;                                                    \
    id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(MPI_COMM_WORLD);                   \
    id.job = JOB_INIT; CFUN(&id);                                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; goto cleanup; }                 \
    if (!verbose) quiet_icntl(id.icntl);                                         \
    if (is_host) { id.n = n; id.nnz = (MUMPS_INT8)nz;                            \
                   id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs; }         \
    id.job = JOB_ANALYZE_FACTOR_SOLVE; CFUN(&id);                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; id.job = JOB_END; CFUN(&id); goto cleanup; } \
    if (is_host) for (int i = 0; i < n; i++) xr[i] = rhs[i].limbs[0];            \
    id.job = JOB_END; CFUN(&id);                                                 \
cleanup:                                                                         \
    free(irn); free(jcn); free(a); free(rhs);                                    \
    return infog1;                                                               \
}

#define GEN_CPLX_DD(NAME, STRUC, CFUN, CPLX)                                      \
static int NAME(const MM *A, const MM *b, int is_host, int sym, int verbose,     \
                double *xr, double *xi) {                                        \
    STRUC id; memset(&id, 0, sizeof id);                                         \
    MUMPS_INT *irn = NULL, *jcn = NULL; CPLX *a = NULL, *rhs = NULL;             \
    int n = 0; long nz = 0; int infog1 = 0;                                      \
    if (is_host) {                                                               \
        n = A->n; nz = A->nnz;                                                   \
        irn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        jcn = malloc(sizeof(MUMPS_INT) * (size_t)nz);                            \
        a   = malloc(sizeof(CPLX) * (size_t)nz);                                 \
        rhs = malloc(sizeof(CPLX) * (size_t)n);                                  \
        for (long k = 0; k < nz; k++) { irn[k] = A->irn[k]; jcn[k] = A->jcn[k];  \
            a[k].r.limbs[0] = A->val[k];                     a[k].r.limbs[1] = 0.0;  \
            a[k].i.limbs[0] = (A->ival ? A->ival[k] : 0.0);  a[k].i.limbs[1] = 0.0; } \
        for (int i = 0; i < n; i++) {                                            \
            rhs[i].r.limbs[0] = b->val[i];                    rhs[i].r.limbs[1] = 0.0;  \
            rhs[i].i.limbs[0] = (b->ival ? b->ival[i] : 0.0); rhs[i].i.limbs[1] = 0.0; } \
    }                                                                            \
    id.par = 1; id.sym = sym;                                                    \
    id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(MPI_COMM_WORLD);                   \
    id.job = JOB_INIT; CFUN(&id);                                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; goto cleanup; }                 \
    if (!verbose) quiet_icntl(id.icntl);                                         \
    if (is_host) { id.n = n; id.nnz = (MUMPS_INT8)nz;                            \
                   id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs; }         \
    id.job = JOB_ANALYZE_FACTOR_SOLVE; CFUN(&id);                                \
    if (id.infog[0] < 0) { infog1 = id.infog[0]; id.job = JOB_END; CFUN(&id); goto cleanup; } \
    if (is_host) for (int i = 0; i < n; i++) { xr[i] = rhs[i].r.limbs[0]; xi[i] = rhs[i].i.limbs[0]; } \
    id.job = JOB_END; CFUN(&id);                                                 \
cleanup:                                                                         \
    free(irn); free(jcn); free(a); free(rhs);                                    \
    return infog1;                                                               \
}

GEN_REAL(solve_s, SMUMPS_STRUC_C, smumps_c, float)
GEN_REAL(solve_d, DMUMPS_STRUC_C, dmumps_c, double)
GEN_CPLX(solve_c, CMUMPS_STRUC_C, cmumps_c, mumps_complex,        float)
GEN_CPLX(solve_z, ZMUMPS_STRUC_C, zmumps_c, mumps_double_complex, double)

/* Extended kind10 (long double) — same macros, wider scalar type. */
GEN_REAL(solve_e, EMUMPS_STRUC_C, emumps_c, long double)
GEN_CPLX(solve_y, YMUMPS_STRUC_C, ymumps_c, mumps_long_double_complex, long double)

/* Extended kind16 (__float128) — scalar cast still works for the value type. */
GEN_REAL(solve_q, QMUMPS_STRUC_C, qmumps_c, mumps_float128)
GEN_CPLX(solve_x, XMUMPS_STRUC_C, xmumps_c, mumps_float128_complex, mumps_float128)

/* Extended multifloats (double-double) — struct value type, dedicated macros. */
GEN_REAL_DD(solve_m, MMUMPS_STRUC_C, mmumps_c, mumps_float64x2)
GEN_CPLX_DD(solve_w, WMUMPS_STRUC_C, wmumps_c, mumps_complex64x2)

static void usage(const char *prog) {
    fprintf(stderr,
        "usage: %s -t <s|c|d|z|e|y|q|x|m|w> [-v] <matrix.mtx> <rhs.mtx> <solution.mtx>\n"
        "  -t  arithmetic (real / complex):\n"
        "        s/c = single        d/z = double        e/y = long double (kind10)\n"
        "        q/x = __float128 (kind16)               m/w = double-double (multifloats)\n"
        "  -v  leave MUMPS diagnostics on (default: silent)\n", prog);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    /* Register the extended-precision custom MPI reductions (idempotent).
     * Needed by q/x and m/w; harmless for the other arithmetics. Done here
     * because MKL's blacs_pinfo_ wins the link and won't register them. */
    quad_mpi_init();
    multifloats_mpi_init();
    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int is_host = (rank == 0);

    /* ── argument parse (every rank sees the same argv) ───────────────*/
    char type = 0; int verbose = 0;
    const char *paths[3] = {0}; int np = 0;
    for (int i = 1; i < argc; i++) {
        if      (strcmp(argv[i], "-t") == 0 && i + 1 < argc) type = argv[++i][0];
        else if (strcmp(argv[i], "-v") == 0)                 verbose = 1;
        else if (np < 3)                                     paths[np++] = argv[i];
        else {
            if (is_host) usage(argv[0]);
            MPI_Finalize(); return 2;
        }
    }
    if (np != 3 || !strchr("scdzeyqxmw", type) || type == 0) {
        if (is_host) usage(argv[0]);
        MPI_Finalize(); return 2;
    }
    const char *mpath = paths[0], *bpath = paths[1], *xpath = paths[2];
    int is_cplx = (strchr("czyxw", type) != NULL);

    /* ── host reads + validates; broadcast an error flag + sym ────────*/
    MM A, b; memset(&A, 0, sizeof A); memset(&b, 0, sizeof b);
    int err = 0, sym = 0;
    if (is_host) {
        if (mm_read(mpath, &A) || mm_read(bpath, &b)) {
            err = 1;
        } else if (!A.is_coordinate) {
            fprintf(stderr, "mmsolve: matrix must be coordinate (sparse)\n"); err = 1;
        } else if (A.m != A.n) {
            fprintf(stderr, "mmsolve: matrix must be square (%d x %d)\n", A.m, A.n); err = 1;
        } else if (b.is_coordinate || b.n != 1 || b.m != A.n) {
            fprintf(stderr, "mmsolve: rhs must be a dense %d x 1 array\n", A.n); err = 1;
        } else if (!is_cplx && (A.is_complex || b.is_complex)) {
            fprintf(stderr, "mmsolve: type '%c' is real but the data is complex\n", type); err = 1;
        } else {
            sym = A.is_symmetric ? 2 : 0;
        }
    }
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (err) { mm_free(&A); mm_free(&b); MPI_Finalize(); return 1; }
    MPI_Bcast(&sym, 1, MPI_INT, 0, MPI_COMM_WORLD);

    /* ── solve ────────────────────────────────────────────────────────*/
    int n = is_host ? A.n : 0;
    double *xr = NULL, *xi = NULL;
    if (is_host) { xr = calloc((size_t)n, sizeof(double));
                   if (is_cplx) xi = calloc((size_t)n, sizeof(double)); }

    int infog = 0;
    switch (type) {
        case 's': infog = solve_s(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'd': infog = solve_d(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'c': infog = solve_c(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'z': infog = solve_z(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'e': infog = solve_e(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'y': infog = solve_y(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'q': infog = solve_q(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'x': infog = solve_x(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'm': infog = solve_m(&A, &b, is_host, sym, verbose, xr, xi); break;
        case 'w': infog = solve_w(&A, &b, is_host, sym, verbose, xr, xi); break;
    }

    int rc = 0;
    if (infog < 0) {
        if (is_host) fprintf(stderr, "mmsolve: MUMPS failed, INFOG(1)=%d\n", infog);
        rc = 1;
    } else if (is_host) {
        if (mm_write_vector(xpath, n, xr, is_cplx ? xi : NULL)) rc = 1;
        else printf("mmsolve: type=%c  n=%d  nnz=%ld  sym=%d  ->  %s\n",
                    type, n, A.nnz, sym, xpath);
    }

    free(xr); free(xi); mm_free(&A); mm_free(&b);
    MPI_Finalize();
    return rc;
}
