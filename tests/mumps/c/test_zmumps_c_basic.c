/*
 * Basic smoke test for the migrated MUMPS C bridge (complex path) —
 * mirror of test_dmumps_c_basic.c with complex (r/i) entries.
 *
 * mumps_double_complex's layout matches the migrated Fortran COMPLEX
 * type's ABI: kind16 → two consecutive __float128 fields, kind10 → two
 * long double fields.  A struct of two scalars with no extra padding
 * is binary-compatible with the corresponding gfortran COMPLEX(KIND).
 *
 * Target portability via test_real_compat.h + TARGET_REAL_*.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "test_real_compat.h"
#include TARGET_COMPLEX_HEADER

#define MUMPS_C       TARGET_COMPLEX_MUMPS_C
#define MUMPS_STRUC_C TARGET_COMPLEX_STRUC_C

#ifdef TEST_TARGET_MULTIFLOATS
extern void multifloats_mpi_init(void);
#endif

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
    printf("  test_zmumps_c_basic [%s] max_rel_err=%s  %s\n",
           case_label, relbuf, passed ? "PASS" : "FAIL");
}

static void report_finalize_c(void)
{
    fprintf(gJson, "  ]\n}\n");
    fclose(gJson);
}

static int report_status_c(void) { return gAnyFail; }

/* Multiplication of test_complex = {r, i}. */
static test_complex cmul(test_complex a, test_complex b)
{
    test_complex c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.r * b.i + a.i * b.r;
    return c;
}

static test_real cabs_real(test_complex z)
{
    return TR_SQRT(z.r * z.r + z.i * z.i);
}

int main(int argc, char **argv)
{
    MUMPS_STRUC_C id = {0};   /* L-5 — see test_dmumps_c_basic.c */
    enum { N = 4 };

    /* Diagonally-dominant complex unsymmetric A. Hand-built so the
     * test is reproducible and the residual stays bounded. */
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

    MUMPS_INT     irn[N*N], jcn[N*N];
    test_complex  a_vals[N*N], rhs[N];

    MPI_Init(&argc, &argv);
#ifdef TEST_TARGET_MULTIFLOATS
    multifloats_mpi_init();
#endif
    {
        int i, j, k = 0;
        for (j = 0; j < N; j++)
            for (i = 0; i < N; i++, k++) {
                irn[k] = i+1; jcn[k] = j+1; a_vals[k] = A[i][j];
            }
        /* b = A * x_true */
        for (i = 0; i < N; i++) {
            test_complex s = { TR_LIT(0.0), TR_LIT(0.0) };
            for (j = 0; j < N; j++) {
                test_complex p = cmul(A[i][j], x_true[j]);
                s.r += p.r; s.i += p.i;
            }
            rhs[i] = s;
        }
    }

    report_init_c("test_zmumps_c_basic", TEST_TARGET_NAME);

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    MUMPS_C(&id);
    if (id.infog[0] < 0) { fprintf(stderr, "JOB=-1 failed: %d\n", id.infog[0]); return 1; }

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    id.n   = N;
    id.nnz = (MUMPS_INT8) (N * N);
    id.irn = irn;
    id.jcn = jcn;
#ifdef TEST_TARGET_MULTIFLOATS
    /* Bridge expects mumps_complex64x2*; widen test_complex (plain
     * doubles) at the boundary and narrow back after the solve. */
    mumps_complex64x2 a_bridge[N*N], rhs_bridge[N];
    for (int k = 0; k < N*N; k++) a_bridge[k]   = tc_widen(a_vals[k]);
    for (int k = 0; k < N;   k++) rhs_bridge[k] = tc_widen(rhs[k]);
    id.a   = a_bridge;
    id.rhs = rhs_bridge;
#else
    id.a   = a_vals;
    id.rhs = rhs;
#endif

    id.job = 6;
    MUMPS_C(&id);
    if (id.infog[0] < 0) {
        fprintf(stderr, "JOB=6 failed: infog=%d %d\n", id.infog[0], id.infog[1]);
        return 1;
    }

#ifdef TEST_TARGET_MULTIFLOATS
    for (int k = 0; k < N; k++) rhs[k] = tc_narrow(rhs_bridge[k]);
#endif

    {
        test_real max_rel = TR_LIT(0.0), denom = TR_LIT(0.0);
        int i;
        for (i = 0; i < N; i++) {
            test_real a = cabs_real(x_true[i]);
            if (a > denom) denom = a;
        }
        for (i = 0; i < N; i++) {
            test_complex d = { rhs[i].r - x_true[i].r, rhs[i].i - x_true[i].i };
            test_real m = cabs_real(d);
            if (m > max_rel) max_rel = m;
        }
        if (denom > TR_MIN) max_rel /= denom;
        report_case_c("n=4", max_rel,
                      TR_LIT(16.0) * (test_real)(N*N*N) * TR_EPS);
    }

    id.job = -2;
    MUMPS_C(&id);

    report_finalize_c();
    MPI_Finalize();
    return report_status_c() ? 1 : 0;
}
