/*
 * Basic smoke test for the migrated MUMPS C bridge (real path).  Mirrors
 * tests/mumps/fortran/test_dmumps_basic.f90 in shape but operates
 * end-to-end through the C interface.
 *
 * Hand-coded (no quad-precision random helpers shared from Fortran)
 * with a small fixed problem: n=4 diagonally-dominant unsymmetric
 * matrix, x_true known, b = A*x_true.  Pass/fail decided by JSON
 * report comparing |x_solve − x_true| / |x_true| against an
 * O(n^3) tolerance scaled by target eps.
 *
 * Target portability is provided by ``test_real_compat.h`` (test_real
 * + TR_LIT/TR_FABS/TR_SQRT/TR_EPS/TR_MIN macros) and the
 * TARGET_REAL_HEADER / TARGET_REAL_MUMPS_C / TARGET_REAL_STRUC_C
 * macros wired in by the cmake glue.
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

#ifdef TEST_TARGET_MULTIFLOATS
extern void multifloats_mpi_init(void);
#endif

/* Same JSON report layout as prec_report.f90 produces, hand-rolled
 * here so the C test side doesn't need to call back into Fortran. */
static FILE *gJson = NULL;
static int gAnyFail = 0;
static int gCaseCount = 0;

static void report_init_c(const char *test_name, const char *target_name)
{
    char filename[256];
    snprintf(filename, sizeof filename, "%s.%s.json", test_name, target_name);
    gJson = fopen(filename, "w");
    if (!gJson) { fprintf(stderr, "report_init_c: cannot open %s\n", filename); exit(1); }
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
    /* Match prec_report.f90: digits = -log10(max_rel), 99.0 sentinel
     * when max_rel == 0. Cast to double for log10 — display only. */
    double digits = (max_rel > (test_real)0)
                        ? -log10((double)max_rel)
                        : 99.0;
    snprintf(digbuf, sizeof digbuf, "%.2f", digits);
    if (gCaseCount > 0) fprintf(gJson, "    ,\n");
    gCaseCount++;
    fprintf(gJson,
            "    {\n"
            "      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %s,\n"
            "      \"tolerance\":   %s,\n"
            "      \"digits\":      %s,\n"
            "      \"passed\":      %s\n"
            "    }\n",
            case_label, relbuf, tolbuf, digbuf,
            passed ? "true" : "false");
    printf("  %s [%s] max_rel_err=%s  %s\n", "test_dmumps_c_basic",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
    gJson = NULL;
}

static int report_status_c(void) { return gAnyFail; }

int main(int argc, char **argv)
{
    /* L-5: zero-init avoids uninit-load of metis_options[40] on the
     * first JOB=-1 call (mumps_c.c reads &id.metis_options before
     * MUMPS_INI_DRIVER overwrites it with defaults). Functionally
     * harmless today but flagged by sanitizers. */
    MUMPS_STRUC_C id = {0};
    int n = 4;
    /* Diagonally-dominant unsymmetric A (column-major triplet). */
    MUMPS_INT  irn[16];
    MUMPS_INT  jcn[16];
    test_real  a_vals[16];
    test_real  x_true[4] = { TR_LIT(1.0), TR_LIT(-2.0), TR_LIT(3.0), TR_LIT(-4.0) };
    test_real  rhs[4];

    MPI_Init(&argc, &argv);
#ifdef TEST_TARGET_MULTIFLOATS
    multifloats_mpi_init();
#endif

    /* Fill A (4x4): diagonal entries large, off-diagonals small.
     * Layout matches dense_to_triplet (column-major flatten). */
    {
        test_real A[4][4] = {
            { TR_LIT( 5.0),  TR_LIT( 0.5),  TR_LIT(-0.25), TR_LIT( 0.1) },
            { TR_LIT( 0.3),  TR_LIT(-6.0),  TR_LIT( 0.5),  TR_LIT( 0.2) },
            { TR_LIT(-0.4),  TR_LIT( 0.2),  TR_LIT( 7.0),  TR_LIT(-0.3) },
            { TR_LIT( 0.1),  TR_LIT(-0.3),  TR_LIT( 0.4),  TR_LIT(-8.0) },
        };
        int i, j, k = 0;
        for (j = 0; j < n; j++) {
            for (i = 0; i < n; i++, k++) {
                irn[k]    = i + 1;
                jcn[k]    = j + 1;
                a_vals[k] = A[i][j];
            }
        }
        /* b = A * x_true */
        for (i = 0; i < n; i++) {
            test_real s = TR_LIT(0.0);
            for (j = 0; j < n; j++) s += A[i][j] * x_true[j];
            rhs[i] = s;
        }
    }

    report_init_c("test_dmumps_c_basic", TEST_TARGET_NAME);

    /* ── Init ─────────────────────────────────────────────────────── */
    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    MUMPS_C(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=-1 failed: infog[0]=%d\n", id.infog[0]);
        return 1;
    }

    /* Silence diagnostic output. */
    id.icntl[0] = -1;  /* error stream */
    id.icntl[1] = -1;  /* diagnostic */
    id.icntl[2] = -1;  /* global info */
    id.icntl[3] =  0;  /* print level */

    /* ── Problem data ────────────────────────────────────────────── */
    id.n    = n;
    id.nnz  = (MUMPS_INT8) (n * n);
    id.irn  = irn;
    id.jcn  = jcn;
#ifdef TEST_TARGET_MULTIFLOATS
    /* Bridge expects mumps_float64x2*; the test data lives in plain
     * double.  Widen at the boundary, narrow back after the solve. */
    mumps_float64x2 a_bridge[16], rhs_bridge[4];
    for (int i = 0; i < n * n; i++) a_bridge[i]   = tr_widen(a_vals[i]);
    for (int i = 0; i < n;     i++) rhs_bridge[i] = tr_widen(rhs[i]);
    id.a    = a_bridge;
    id.rhs  = rhs_bridge;
#else
    id.a    = a_vals;
    id.rhs  = rhs;
#endif

    /* ── Solve ───────────────────────────────────────────────────── */
    id.job = 6;
    MUMPS_C(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=6 failed: infog[0]=%d, infog[1]=%d\n",
                id.infog[0], id.infog[1]);
        return 1;
    }

#ifdef TEST_TARGET_MULTIFLOATS
    for (int i = 0; i < n; i++) rhs[i] = tr_narrow(rhs_bridge[i]);
#endif

    /* On exit, rhs[] holds the solution. */
    {
        test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
        int i;
        for (i = 0; i < n; i++) {
            test_real ax = TR_FABS(x_true[i]);
            if (ax > denom) denom = ax;
        }
        for (i = 0; i < n; i++) {
            test_real d = TR_FABS(rhs[i] - x_true[i]);
            if (d > max_rel) max_rel = d;
        }
        if (denom > TR_MIN) max_rel /= denom;

        test_real tol = TR_LIT(16.0) * (test_real) (n * n * n) * TR_EPS;
        report_case_c("n=4", max_rel, tol);
    }

    /* ── Cleanup ─────────────────────────────────────────────────── */
    id.job = -2;
    MUMPS_C(&id);

    report_finalize_c();
    MPI_Finalize();
    return report_status_c() ? 1 : 0;
}
