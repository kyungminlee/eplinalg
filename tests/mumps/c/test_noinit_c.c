/*
 * Sequential (libmpiseq) no-initialization solve check — C bridge side.
 *
 * Mirror of tests/mumps/fortran/test_noinit_f.f90, driving the migrated
 * MUMPS C bridge for BOTH supported types of this target (real +
 * complex) with NEITHER MPI_Init NOR the per-target MPI-handle
 * registration (test_target_mpi_init / multifloats_mpi_init /
 * quad_mpi_init). This is exactly how a plain sequential C consumer uses
 * the release: link libmpiseq, call ?mumps_c, never touch MPI.
 *
 * On multifloats the MPI_FLOAT64X2 / MPI_COMPLEX64X2 datatype handles now
 * default to the libmpiseq derived-type sentinel, so libseq's MUMPS_COPY
 * dispatches correctly without multifloats_mpi_init having run (see
 * runtime/multifloats-mpi/multifloats_mpi.cpp). kind16 / kind10 never
 * needed init for a single-rank solve.
 *
 * Built ONLY in the _seq variant (see tests/mumps/CMakeLists.txt): the
 * real-MPI variant legitimately requires MPI_Init.
 *
 * Each type solves a small n=4 dense system A x = b (JOB=-1 -> 6 -> -2)
 * and checks BOTH the solution accuracy and the residual norm against an
 * O(n^3) eps tolerance. Deliberately no <mpi.h>: id.comm_fortran = 0 is
 * ignored by libmpiseq, proving the path needs nothing from MPI.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_real_compat.h"   /* pulls both TARGET headers + types/macros */

#define REAL_MUMPS_C     TARGET_REAL_MUMPS_C
#define REAL_STRUC_C     TARGET_REAL_STRUC_C
#define COMPLEX_MUMPS_C  TARGET_COMPLEX_MUMPS_C
#define COMPLEX_STRUC_C  TARGET_COMPLEX_STRUC_C

enum { N = 4 };

/* ── JSON report (shared across the four cases this program emits) ──── */
static FILE *gJson = NULL;
static int gAnyFail = 0, gCaseCount = 0;

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
    double digits = (max_rel > (test_real)0) ? -log10((double)max_rel) : 99.0;
    snprintf(digbuf, sizeof digbuf, "%.2f", digits);
    if (gCaseCount > 0) fprintf(gJson, "    ,\n");
    gCaseCount++;
    fprintf(gJson,
            "    {\n      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %s,\n      \"tolerance\":   %s,\n"
            "      \"digits\":      %s,\n      \"passed\":      %s\n    }\n",
            case_label, relbuf, tolbuf, digbuf, passed ? "true" : "false");
    printf("  test_noinit_c [%s] max_rel_err=%s  %s\n",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void) { fprintf(gJson, "  ]\n}\n"); fclose(gJson); }
static int  report_status_c(void)   { return gAnyFail; }

/* ── complex scalar helpers on test_complex = {r, i} ───────────────── */
static test_complex cmul(test_complex a, test_complex b)
{
    test_complex c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.r * b.i + a.i * b.r;
    return c;
}
static test_real cabs_real(test_complex z) { return TR_SQRT(z.r * z.r + z.i * z.i); }

/* ── real path ─────────────────────────────────────────────────────── */
static void solve_real(void)
{
    REAL_STRUC_C id = {0};   /* L-5 zero-init — see test_dmumps_c_basic.c */
    test_real A[N][N] = {
        { TR_LIT( 5.0),  TR_LIT( 0.5),  TR_LIT(-0.25), TR_LIT( 0.1) },
        { TR_LIT( 0.3),  TR_LIT(-6.0),  TR_LIT( 0.5),  TR_LIT( 0.2) },
        { TR_LIT(-0.4),  TR_LIT( 0.2),  TR_LIT( 7.0),  TR_LIT(-0.3) },
        { TR_LIT( 0.1),  TR_LIT(-0.3),  TR_LIT( 0.4),  TR_LIT(-8.0) },
    };
    test_real x_true[N] = { TR_LIT(1.0), TR_LIT(-2.0), TR_LIT(3.0), TR_LIT(-4.0) };
    MUMPS_INT irn[N*N], jcn[N*N];
    test_real a_vals[N*N], bvec[N], rhs[N];
    int i, j, k = 0;

    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++, k++) {
            irn[k] = i + 1; jcn[k] = j + 1; a_vals[k] = A[i][j];
        }
    for (i = 0; i < N; i++) {
        test_real s = TR_LIT(0.0);
        for (j = 0; j < N; j++) s += A[i][j] * x_true[j];
        bvec[i] = s; rhs[i] = s;
    }

    /* NO MPI_Init, NO test_target_mpi_init. comm_fortran ignored by libmpiseq. */
    id.par = 1; id.sym = 0; id.comm_fortran = 0; id.job = -1;
    REAL_MUMPS_C(&id);
    if (id.infog[0] < 0) { fprintf(stderr, "real JOB=-1 failed: %d\n", id.infog[0]); exit(1); }
    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

#ifdef TEST_TARGET_MULTIFLOATS
    mumps_float64x2 a_bridge[N*N], rhs_bridge[N];
    for (k = 0; k < N*N; k++) a_bridge[k]   = tr_widen(a_vals[k]);
    for (k = 0; k < N;   k++) rhs_bridge[k] = tr_widen(rhs[k]);
    id.a = a_bridge; id.rhs = rhs_bridge;
#else
    id.a = a_vals; id.rhs = rhs;
#endif
    id.n = N; id.nnz = (MUMPS_INT8)(N * N); id.irn = irn; id.jcn = jcn;

    id.job = 6;
    REAL_MUMPS_C(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "real JOB=6 failed: %d %d\n", id.infog[0], id.infog[1]); exit(1);
    }
#ifdef TEST_TARGET_MULTIFLOATS
    for (k = 0; k < N; k++) rhs[k] = tr_narrow(rhs_bridge[k]);
#endif

    /* solution accuracy */
    {
        test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
        for (i = 0; i < N; i++) { test_real a = TR_FABS(x_true[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) { test_real d = TR_FABS(rhs[i] - x_true[i]); if (d > max_rel) max_rel = d; }
        if (denom > TR_MIN) max_rel /= denom;
        report_case_c("real solution", max_rel, TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
    }
    /* residual  max|b - A x_solve| / max|b| */
    {
        test_real max_res = TR_LIT(0.0), denom = TR_LIT(0.0);
        for (i = 0; i < N; i++) { test_real a = TR_FABS(bvec[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) {
            test_real s = TR_LIT(0.0);
            for (j = 0; j < N; j++) s += A[i][j] * rhs[j];
            test_real d = TR_FABS(bvec[i] - s);
            if (d > max_res) max_res = d;
        }
        if (denom > TR_MIN) max_res /= denom;
        report_case_c("real residual", max_res, TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
    }

    id.job = -2;
    REAL_MUMPS_C(&id);
}

/* ── complex path ──────────────────────────────────────────────────── */
static void solve_complex(void)
{
    COMPLEX_STRUC_C id = {0};
    test_complex A[N][N] = {
        {{TR_LIT( 5.0), TR_LIT( 0.5)}, {TR_LIT( 0.5), TR_LIT( 0.1)},
         {TR_LIT(-0.25),TR_LIT( 0.05)},{TR_LIT( 0.1), TR_LIT( 0.0)}},
        {{TR_LIT( 0.3), TR_LIT( 0.0)}, {TR_LIT(-6.0), TR_LIT(-0.4)},
         {TR_LIT( 0.5), TR_LIT( 0.1)}, {TR_LIT( 0.2), TR_LIT(-0.1)}},
        {{TR_LIT(-0.4), TR_LIT( 0.1)}, {TR_LIT( 0.2), TR_LIT(-0.05)},
         {TR_LIT( 7.0), TR_LIT( 0.3)}, {TR_LIT(-0.3), TR_LIT( 0.2)}},
        {{TR_LIT( 0.1), TR_LIT(-0.05)},{TR_LIT(-0.3), TR_LIT( 0.1)},
         {TR_LIT( 0.4), TR_LIT( 0.0)}, {TR_LIT(-8.0), TR_LIT( 0.5)}},
    };
    test_complex x_true[N] = {
        {TR_LIT( 1.0), TR_LIT(-0.5)}, {TR_LIT(-2.0), TR_LIT( 1.0)},
        {TR_LIT( 3.0), TR_LIT(-1.5)}, {TR_LIT(-4.0), TR_LIT( 2.0)},
    };
    MUMPS_INT irn[N*N], jcn[N*N];
    test_complex a_vals[N*N], bvec[N], rhs[N];
    int i, j, k = 0;

    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++, k++) {
            irn[k] = i + 1; jcn[k] = j + 1; a_vals[k] = A[i][j];
        }
    for (i = 0; i < N; i++) {
        test_complex s = { TR_LIT(0.0), TR_LIT(0.0) };
        for (j = 0; j < N; j++) { test_complex p = cmul(A[i][j], x_true[j]); s.r += p.r; s.i += p.i; }
        bvec[i] = s; rhs[i] = s;
    }

    id.par = 1; id.sym = 0; id.comm_fortran = 0; id.job = -1;
    COMPLEX_MUMPS_C(&id);
    if (id.infog[0] < 0) { fprintf(stderr, "complex JOB=-1 failed: %d\n", id.infog[0]); exit(1); }
    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

#ifdef TEST_TARGET_MULTIFLOATS
    mumps_complex64x2 a_bridge[N*N], rhs_bridge[N];
    for (k = 0; k < N*N; k++) a_bridge[k]   = tc_widen(a_vals[k]);
    for (k = 0; k < N;   k++) rhs_bridge[k] = tc_widen(rhs[k]);
    id.a = a_bridge; id.rhs = rhs_bridge;
#else
    id.a = a_vals; id.rhs = rhs;
#endif
    id.n = N; id.nnz = (MUMPS_INT8)(N * N); id.irn = irn; id.jcn = jcn;

    id.job = 6;
    COMPLEX_MUMPS_C(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "complex JOB=6 failed: %d %d\n", id.infog[0], id.infog[1]); exit(1);
    }
#ifdef TEST_TARGET_MULTIFLOATS
    for (k = 0; k < N; k++) rhs[k] = tc_narrow(rhs_bridge[k]);
#endif

    {
        test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
        for (i = 0; i < N; i++) { test_real a = cabs_real(x_true[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) {
            test_complex d = { rhs[i].r - x_true[i].r, rhs[i].i - x_true[i].i };
            test_real m = cabs_real(d); if (m > max_rel) max_rel = m;
        }
        if (denom > TR_MIN) max_rel /= denom;
        report_case_c("complex solution", max_rel, TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
    }
    {
        test_real max_res = TR_LIT(0.0), denom = TR_LIT(0.0);
        for (i = 0; i < N; i++) { test_real a = cabs_real(bvec[i]); if (a > denom) denom = a; }
        for (i = 0; i < N; i++) {
            test_complex s = { TR_LIT(0.0), TR_LIT(0.0) };
            for (j = 0; j < N; j++) { test_complex p = cmul(A[i][j], rhs[j]); s.r += p.r; s.i += p.i; }
            test_complex d = { bvec[i].r - s.r, bvec[i].i - s.i };
            test_real m = cabs_real(d); if (m > max_res) max_res = m;
        }
        if (denom > TR_MIN) max_res /= denom;
        report_case_c("complex residual", max_res, TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
    }

    id.job = -2;
    COMPLEX_MUMPS_C(&id);
}

int main(void)
{
    report_init_c("test_noinit_c", TEST_TARGET_NAME);
    solve_real();
    solve_complex();
    report_finalize_c();
    return report_status_c() ? 1 : 0;
}
