/*
 * SYM coverage via the migrated MUMPS C bridge (real path) — mirrors
 * tests/mumps/fortran/test_dmumps_sym.f90 with hand-coded n=4
 * problems for each SYM ∈ {0, 1, 2}.
 *
 * SYM=0  :  general unsymmetric A (LU)
 * SYM=1  :  symmetric positive-definite A = X*X^T + n*I (Cholesky)
 * SYM=2  :  general symmetric A = (R + R^T)/2 + n*I (Bunch-Kaufman)
 *
 * For SYM>0 only the upper triangle is sent in the triplet form; the
 * solver infers the lower triangle.
 *
 * Target portability via test_real_compat.h + TARGET_REAL_*.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "test_real_compat.h"
#include TARGET_REAL_HEADER

#define MUMPS_C       TARGET_REAL_MUMPS_C
#define MUMPS_STRUC_C TARGET_REAL_STRUC_C


static FILE *gJson = NULL;
static int gAnyFail = 0;
static int gCaseCount = 0;

static void report_init_c(const char *test_name, const char *target_name)
{
    char filename[256];
    snprintf(filename, sizeof filename, "%s.%s.json", test_name, target_name);
    gJson = fopen(filename, "w");
    if (!gJson) { fprintf(stderr, "report_init: cannot open %s\n", filename); exit(1); }
    fprintf(gJson, "{\n  \"routine\": \"%s\",\n  \"target\":  \"%s\",\n  \"cases\": [\n",
            test_name, target_name);
}

static void report_case_c(const char *case_label, test_real max_rel, test_real tol)
{
    char relbuf[64], tolbuf[64], digbuf[64];
    int passed = (max_rel <= tol);
    if (!passed) gAnyFail = 1;
    test_real_snprintf(relbuf, sizeof relbuf, max_rel);
    test_real_snprintf(tolbuf, sizeof tolbuf, tol);
    double digits = (max_rel > (test_real)0)
                        ? -log10((double)max_rel)
                        : 99.0;
    snprintf(digbuf, sizeof digbuf, "%.2f", digits);
    if (gCaseCount > 0) fprintf(gJson, "    ,\n");
    gCaseCount++;
    fprintf(gJson,
            "    {\n      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %s,\n      \"tolerance\":   %s,\n"
            "      \"digits\":      %s,\n"
            "      \"passed\":      %s\n    }\n",
            case_label, relbuf, tolbuf, digbuf,
            passed ? "true" : "false");
    printf("  test_dmumps_c_sym [%s] max_rel_err=%s  %s\n",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
}

static int report_status_c(void) { return gAnyFail; }

/* Solve A*x = b via the migrated MUMPS bridge, given a triplet
 * (irn, jcn, a) of nz entries and an n-vector b (overwritten with x
 * on exit). On multifloats the test data is plain double — widen
 * to the bridge's mumps_float64x2 buffer at the boundary, narrow
 * back into rhs after the solve. The compile-time-bounded N*(N+1)/2
 * upper-triangle worst-case is N*N; the bridge buffers are sized
 * accordingly. */
static int mumps_solve(int sym, MUMPS_INT n, MUMPS_INT8 nnz,
                        MUMPS_INT *irn, MUMPS_INT *jcn,
                        test_real *a_vals, test_real *rhs)
{
    MUMPS_STRUC_C id = {0};   /* L-5 — see test_dmumps_c_basic.c */
    /* Centralized input: matrix + RHS on the host only; the centralized
     * solution comes back in rhs[] on the host only. MUMPS calls stay
     * collective. Mirrors test_dmumps_basic.f90. */
    int myid, is_host;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    is_host = (myid == 0);
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    MUMPS_C(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

#ifdef TEST_TARGET_MULTIFLOATS
    /* Declared at function scope so the buffers outlive the host guard. */
    mumps_float64x2 a_bridge[16], rhs_bridge[4];
#endif
    if (is_host) {
        id.n   = n;
        id.nnz = nnz;
        id.irn = irn;
        id.jcn = jcn;
#ifdef TEST_TARGET_MULTIFLOATS
        for (MUMPS_INT8 k = 0; k < nnz; k++) a_bridge[k]   = tr_widen(a_vals[k]);
        for (MUMPS_INT  k = 0; k < n;   k++) rhs_bridge[k] = tr_widen(rhs[k]);
        id.a   = a_bridge;
        id.rhs = rhs_bridge;
#else
        id.a   = a_vals;
        id.rhs = rhs;
#endif
    }
    id.job = 6;
    MUMPS_C(&id);
    int code = id.infog[0];

#ifdef TEST_TARGET_MULTIFLOATS
    if (is_host)
        for (MUMPS_INT k = 0; k < n; k++) rhs[k] = tr_narrow(rhs_bridge[k]);
#endif

    id.job = -2;
    MUMPS_C(&id);
    return code;
}

int main(int argc, char **argv)
{
    enum { N = 4 };
    test_real x_true[N] = { TR_LIT(1.0), TR_LIT(-2.0), TR_LIT(3.0), TR_LIT(-4.0) };
    int myid, is_host;

    MPI_Init(&argc, &argv);
    test_target_mpi_init();
    /* Centralized input/solution: matrix build, solution check and JSON
     * report are host-guarded; mumps_solve() stays collective on every
     * rank. Mirrors test_dmumps_basic.f90. */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    is_host = (myid == 0);
    if (is_host) report_init_c("test_dmumps_c_sym", TEST_TARGET_NAME);

    /* ── SYM=0 — general unsymmetric ─────────────────────────────── */
    {
        test_real A[N][N] = {
            { TR_LIT( 5.0),  TR_LIT( 0.5),  TR_LIT(-0.25), TR_LIT( 0.1) },
            { TR_LIT( 0.3),  TR_LIT(-6.0),  TR_LIT( 0.5),  TR_LIT( 0.2) },
            { TR_LIT(-0.4),  TR_LIT( 0.2),  TR_LIT( 7.0),  TR_LIT(-0.3) },
            { TR_LIT( 0.1),  TR_LIT(-0.3),  TR_LIT( 0.4),  TR_LIT(-8.0) },
        };
        MUMPS_INT  irn[N*N], jcn[N*N];
        test_real  vals[N*N], rhs[N];
        int i, j, k = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++, k++) {
                irn[k] = i+1; jcn[k] = j+1; vals[k] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            test_real s = TR_LIT(0.0);
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = mumps_solve(0, N, (MUMPS_INT8)(N*N), irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=0 failed: %d\n", code); return 1; }
        if (is_host) {
            test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
            for (i = 0; i < N; i++) { test_real a = TR_FABS(x_true[i]); if (a > denom) denom = a; }
            for (i = 0; i < N; i++) { test_real d = TR_FABS(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
            if (denom > TR_MIN) max_rel /= denom;
            report_case_c("sym=0", max_rel,
                          TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
        }
    }

    /* ── SYM=1 — symmetric positive-definite ─────────────────────── */
    {
        /* A = X*X^T + N*I — guaranteed SPD. Hand-build it. */
        test_real X[N][N] = {
            { TR_LIT( 1.0),  TR_LIT( 0.5),  TR_LIT(-0.25), TR_LIT( 0.1) },
            { TR_LIT( 0.3),  TR_LIT( 0.8),  TR_LIT( 0.5),  TR_LIT( 0.2) },
            { TR_LIT(-0.4),  TR_LIT( 0.2),  TR_LIT( 0.7),  TR_LIT(-0.3) },
            { TR_LIT( 0.1),  TR_LIT(-0.3),  TR_LIT( 0.4),  TR_LIT( 0.6) },
        };
        test_real A[N][N] = { 0 };
        int i, j, k;
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                test_real s = TR_LIT(0.0);
                for (k = 0; k < N; k++) s += X[i][k] * X[j][k];
                A[i][j] = s + (i == j ? (test_real) N : TR_LIT(0.0));
            }
        /* Upper-triangle triplet for SYM=1 */
        MUMPS_INT  irn[N*(N+1)/2], jcn[N*(N+1)/2];
        test_real  vals[N*(N+1)/2], rhs[N];
        int idx = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i <= j; i++, idx++) {
                irn[idx] = i+1; jcn[idx] = j+1; vals[idx] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            test_real s = TR_LIT(0.0);
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = mumps_solve(1, N, (MUMPS_INT8) idx, irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=1 failed: %d\n", code); return 1; }
        if (is_host) {
            test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
            for (i = 0; i < N; i++) { test_real a = TR_FABS(x_true[i]); if (a > denom) denom = a; }
            for (i = 0; i < N; i++) { test_real d = TR_FABS(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
            if (denom > TR_MIN) max_rel /= denom;
            report_case_c("sym=1", max_rel,
                          TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
        }
    }

    /* ── SYM=2 — general symmetric ───────────────────────────────── */
    {
        test_real A[N][N] = {
            { TR_LIT( 6.0),  TR_LIT( 0.5),  TR_LIT(-0.25), TR_LIT( 0.1) },
            { TR_LIT( 0.5),  TR_LIT(-7.0),  TR_LIT( 0.5),  TR_LIT( 0.2) },
            { TR_LIT(-0.25), TR_LIT( 0.5),  TR_LIT( 8.0),  TR_LIT(-0.3) },
            { TR_LIT( 0.1),  TR_LIT( 0.2),  TR_LIT(-0.3),  TR_LIT(-9.0) },
        };  /* hand-symmetric */
        MUMPS_INT  irn[N*(N+1)/2], jcn[N*(N+1)/2];
        test_real  vals[N*(N+1)/2], rhs[N];
        int i, j, idx = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i <= j; i++, idx++) {
                irn[idx] = i+1; jcn[idx] = j+1; vals[idx] = A[i][j];
            }
        for (i = 0; i < N; i++) {
            test_real s = TR_LIT(0.0);
            for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
        int code = mumps_solve(2, N, (MUMPS_INT8) idx, irn, jcn, vals, rhs);
        if (code < 0) { fprintf(stderr, "SYM=2 failed: %d\n", code); return 1; }
        if (is_host) {
            test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
            for (i = 0; i < N; i++) { test_real a = TR_FABS(x_true[i]); if (a > denom) denom = a; }
            for (i = 0; i < N; i++) { test_real d = TR_FABS(rhs[i]-x_true[i]); if (d > max_rel) max_rel = d; }
            if (denom > TR_MIN) max_rel /= denom;
            report_case_c("sym=2", max_rel,
                          TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
        }
    }

    if (is_host) report_finalize_c();
    MPI_Finalize();
    return is_host ? (report_status_c() ? 1 : 0) : 0;
}
