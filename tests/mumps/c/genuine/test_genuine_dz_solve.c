/*
 * Genuine double-precision MUMPS C-bridge solve test.
 *
 * Unlike the ``d``/``z``-named tests under ``c/`` — which the build
 * remaps (LIB_PREFIX / LIB_PREFIX_COMPLEX) onto the MIGRATED bridge and
 * therefore exercise q/x (kind16), e/y (kind10) or m/w (multifloats) at
 * the target's extended precision — this test calls the GENUINE
 * ``dmumps_c`` / ``zmumps_c`` entry points backed by the pristine
 * upstream double-precision libdmumps / libzmumps. It is the only
 * per-stage ctest that verifies the genuine double bridge; residuals
 * come back at plain-``double`` accuracy (~1e-16), not the migrated
 * width.
 *
 * A single executable drives both a real solve (dmumps_c) and a complex
 * solve (zmumps_c): genuine d and z are different arithmetic sharing
 * ``dzmumps_common``, so they coexist in one program (whereas genuine
 * and migrated collide on the unprefixed Fortran commons — see
 * tests/mumps/CMakeLists.txt).
 *
 * Centralized-input MUMPS: matrix + RHS supplied on the host only, the
 * centralized solution returned in id.rhs on the host only. MUMPS calls
 * stay collective on every rank; matrix population, the residual check
 * and the JSON report are host-guarded so the test is correct under
 * mpiexec -n >1 as well as MPISEQ. Mirrors test_dmumps_c_basic.c.
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "dmumps_c.h"
#include "zmumps_c.h"

enum { N = 4 };

static FILE *gJson = NULL;
static int gAnyFail = 0, gCaseCount = 0;

static void report_init(const char *test_name, const char *target_name)
{
    char filename[256];
    snprintf(filename, sizeof filename, "%s.%s.json", test_name, target_name);
    gJson = fopen(filename, "w");
    if (!gJson) { fprintf(stderr, "report_init: cannot open %s\n", filename); exit(1); }
    fprintf(gJson, "{\n  \"routine\": \"%s\",\n  \"target\":  \"%s\",\n  \"cases\": [\n",
            test_name, target_name);
}

static void report_case(const char *case_label, double max_rel, double tol)
{
    int passed = (max_rel <= tol);
    double digits = (max_rel > 0.0) ? -log10(max_rel) : 99.0;
    if (!passed) gAnyFail = 1;
    if (gCaseCount > 0) fprintf(gJson, "    ,\n");
    gCaseCount++;
    fprintf(gJson,
            "    {\n      \"case\":        \"%s\",\n"
            "      \"max_rel_err\": %.6e,\n      \"tolerance\":   %.6e,\n"
            "      \"digits\":      %.2f,\n"
            "      \"passed\":      %s\n    }\n",
            case_label, max_rel, tol, digits, passed ? "true" : "false");
    printf("  test_genuine_dz_solve [%s] max_rel_err=%.6e  %s\n",
           case_label, max_rel, passed ? "PASS" : "FAIL");
}

static void report_finalize(void) { fprintf(gJson, "  ]\n}\n"); fclose(gJson); }

/* Diagonally-dominant 4x4 A (shared real skeleton). */
static const double Ar[N][N] = {
    { 5.0,  0.5, -0.25, 0.1},
    { 0.3, -6.0,  0.5,  0.2},
    {-0.4,  0.2,  7.0, -0.3},
    { 0.1, -0.3,  0.4, -8.0},
};

/* ── real solve via genuine dmumps_c ──────────────────────────────── */
static int solve_real(int is_host, double *rel_out, double *tol_out)
{
    DMUMPS_STRUC_C id;
    memset(&id, 0, sizeof id);

    MUMPS_INT irn[N*N], jcn[N*N];
    double    a[N*N], rhs[N];
    double    xt[N] = {1.0, -2.0, 3.0, -4.0};
    int i, j, k = 0;

    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++, k++) { irn[k] = i+1; jcn[k] = j+1; a[k] = Ar[i][j]; }
    for (i = 0; i < N; i++) { double s = 0.0; for (j = 0; j < N; j++) s += Ar[i][j]*xt[j]; rhs[i] = s; }

    id.par = 1; id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1; dmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];
    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    if (is_host) {
        id.n = N; id.nnz = (MUMPS_INT8)(N*N);
        id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs;
    }
    id.job = 6; dmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];

    if (is_host) {
        double mr = 0.0, den = 0.0;
        for (i = 0; i < N; i++) { double ax = fabs(xt[i]); if (ax > den) den = ax; }
        for (i = 0; i < N; i++) { double d = fabs(rhs[i] - xt[i]); if (d > mr) mr = d; }
        if (den > 0.0) mr /= den;
        *rel_out = mr;
        *tol_out = 16.0 * (double)(N*N*N) * DBL_EPSILON;
    }
    id.job = -2; dmumps_c(&id);
    return 0;
}

/* ── complex solve via genuine zmumps_c ───────────────────────────── */
static int solve_complex(int is_host, double *rel_out, double *tol_out)
{
    ZMUMPS_STRUC_C id;
    memset(&id, 0, sizeof id);

    MUMPS_INT            irn[N*N], jcn[N*N];
    mumps_double_complex a[N*N], rhs[N], xt[N];
    static const double  Ai[N][N] = {
        {0.2, 0.05, 0.0,  0.01},
        {0.03,0.3, -0.05, 0.0 },
        {0.0, 0.02,-0.4,  0.03},
        {0.01,0.0,  0.04, 0.5 },
    };
    mumps_double_complex A[N][N];
    int i, j, k = 0;

    xt[0].r = 1.0;  xt[0].i = 1.0;
    xt[1].r = -2.0; xt[1].i = 0.5;
    xt[2].r = 3.0;  xt[2].i = -1.0;
    xt[3].r = -4.0; xt[3].i = 2.0;
    for (i = 0; i < N; i++) for (j = 0; j < N; j++) { A[i][j].r = Ar[i][j]; A[i][j].i = Ai[i][j]; }
    for (j = 0; j < N; j++)
        for (i = 0; i < N; i++, k++) { irn[k] = i+1; jcn[k] = j+1; a[k] = A[i][j]; }
    for (i = 0; i < N; i++) {
        double sr = 0.0, si = 0.0;
        for (j = 0; j < N; j++) {
            sr += A[i][j].r*xt[j].r - A[i][j].i*xt[j].i;
            si += A[i][j].r*xt[j].i + A[i][j].i*xt[j].r;
        }
        rhs[i].r = sr; rhs[i].i = si;
    }

    id.par = 1; id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1; zmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];
    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    if (is_host) {
        id.n = N; id.nnz = (MUMPS_INT8)(N*N);
        id.irn = irn; id.jcn = jcn; id.a = a; id.rhs = rhs;
    }
    id.job = 6; zmumps_c(&id);
    if (id.infog[0] < 0) return id.infog[0];

    if (is_host) {
        double mr = 0.0, den = 0.0;
        for (i = 0; i < N; i++) { double ax = xt[i].r*xt[i].r + xt[i].i*xt[i].i; if (ax > den) den = ax; }
        for (i = 0; i < N; i++) {
            double dr = rhs[i].r - xt[i].r, di = rhs[i].i - xt[i].i;
            double d = dr*dr + di*di; if (d > mr) mr = d;
        }
        double rel = (den > 0.0) ? mr/den : mr;
        *rel_out = sqrt(rel);
        *tol_out = 16.0 * (double)(N*N*N) * DBL_EPSILON;
    }
    id.job = -2; zmumps_c(&id);
    return 0;
}

int main(int argc, char **argv)
{
    int myid, is_host, code;
    double rel = 0.0, tol = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    is_host = (myid == 0);

    if (is_host) report_init("test_genuine_dz_solve", "genuine");

    code = solve_real(is_host, &rel, &tol);
    if (code < 0) { fprintf(stderr, "genuine dmumps_c failed, infog=%d\n", code); MPI_Finalize(); return 1; }
    if (is_host) report_case("real-n=4", rel, tol);

    code = solve_complex(is_host, &rel, &tol);
    if (code < 0) { fprintf(stderr, "genuine zmumps_c failed, infog=%d\n", code); MPI_Finalize(); return 1; }
    if (is_host) report_case("complex-n=4", rel, tol);

    if (is_host) report_finalize();
    MPI_Finalize();
    return is_host ? (gAnyFail ? 1 : 0) : 0;
}
