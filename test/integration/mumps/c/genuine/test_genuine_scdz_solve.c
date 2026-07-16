/*
 * Genuine all-precision MUMPS C-bridge solve test (s, c, d, z).
 *
 * One executable links every genuine arithmetic library shipped in the
 * release — libsmumps / libcmumps / libdmumps / libzmumps (via their
 * C bridges smumps_c / cmumps_c / dmumps_c / zmumps_c) — and drives a
 * real and a complex solve through each. It is the direct answer to
 * "link all precision libs to a single executable which tests them all":
 * the four genuine arithmetics share the one unpromoted ``mumps_common``
 * (this is what retiring the separate ``dzmumps_common`` unlocked), so
 * they coexist in a single program. Each arithmetic's residual returns
 * at its own native accuracy — single (~1e-6) for s/c, double (~1e-15)
 * for d/z.
 *
 * The migrated extended bridge (e/y, q/x or m/w) is deliberately NOT
 * linked here: its bridge is mumps_c.c compiled at MUMPS_ARITH_d and so
 * emits the same dmumps_* C callback helpers as the genuine d bridge —
 * a C-symbol clash independent of the (now shared) Fortran common. The
 * extended arithmetic is exercised by the remapped c/ tests instead.
 *
 * Centralized-input MUMPS: matrix + RHS supplied on the host only, the
 * centralized solution returned in id.rhs on the host only. MUMPS calls
 * stay collective on every rank; matrix population, the residual check
 * and the JSON report are host-guarded so the test is correct under
 * mpiexec -n >1 as well as MPISEQ. Mirrors test_genuine_dz_solve.c.
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#include "smumps_c.h"
#include "cmumps_c.h"
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
    printf("  test_genuine_scdz_solve [%-14s] max_rel_err=%.6e  %s\n",
           case_label, max_rel, passed ? "PASS" : "FAIL");
}

static void report_finalize(void) { fprintf(gJson, "  ]\n}\n"); fclose(gJson); }

/* Diagonally-dominant 4x4 A (shared real skeleton, exact in binary fp). */
static const double Ar[N][N] = {
    { 5.0,  0.5, -0.25, 0.1},
    { 0.3, -6.0,  0.5,  0.2},
    {-0.4,  0.2,  7.0, -0.3},
    { 0.1, -0.3,  0.4, -8.0},
};
/* Imaginary skeleton for the complex arithmetics. */
static const double Ai[N][N] = {
    {0.2, 0.05, 0.0,  0.01},
    {0.03,0.3, -0.05, 0.0 },
    {0.0, 0.02,-0.4,  0.03},
    {0.01,0.0,  0.04, 0.5 },
};

/* ── real solve, generated per genuine real arithmetic ────────────────
 * REAL is the working type (float for s, double for d); STRUC / CFUN /
 * EPS pick the arithmetic. The residual is measured and reduced in
 * double regardless, so the tolerance is the only precision-specific
 * knob. */
#define GEN_REAL_SOLVE(NAME, STRUC, CFUN, REAL, EPS)                        \
static int NAME(int is_host, double *rel_out, double *tol_out)              \
{                                                                          \
    STRUC id; memset(&id, 0, sizeof id);                                   \
    MUMPS_INT irn[N*N], jcn[N*N];                                          \
    REAL a[N*N], rhs[N];                                                   \
    REAL xt[N] = {(REAL)1.0, (REAL)-2.0, (REAL)3.0, (REAL)-4.0};           \
    int i, j, k = 0;                                                       \
    for (j = 0; j < N; j++)                                                \
        for (i = 0; i < N; i++, k++) { irn[k]=i+1; jcn[k]=j+1; a[k]=(REAL)Ar[i][j]; } \
    for (i = 0; i < N; i++) { REAL s=(REAL)0; for (j=0;j<N;j++) s+=(REAL)Ar[i][j]*xt[j]; rhs[i]=s; } \
    id.par = 1; id.sym = 0;                                                \
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);                        \
    id.job = -1; CFUN(&id);                                                \
    if (id.infog[0] < 0) return id.infog[0];                              \
    id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;        \
    if (is_host) { id.n=N; id.nnz=(MUMPS_INT8)(N*N); id.irn=irn; id.jcn=jcn; id.a=a; id.rhs=rhs; } \
    id.job = 6; CFUN(&id);                                                 \
    if (id.infog[0] < 0) return id.infog[0];                              \
    if (is_host) {                                                         \
        double mr=0.0, den=0.0;                                            \
        for (i=0;i<N;i++){ double ax=fabs((double)xt[i]); if(ax>den) den=ax; } \
        for (i=0;i<N;i++){ double d=fabs((double)rhs[i]-(double)xt[i]); if(d>mr) mr=d; } \
        if (den>0.0) mr/=den;                                              \
        *rel_out=mr; *tol_out=16.0*(double)(N*N*N)*(EPS);                  \
    }                                                                      \
    id.job = -2; CFUN(&id);                                                \
    return 0;                                                              \
}

/* ── complex solve, generated per genuine complex arithmetic ──────────
 * CPLX is the interleaved {r,i} struct (mumps_complex for c,
 * mumps_double_complex for z); REAL is its component type. */
#define GEN_CPLX_SOLVE(NAME, STRUC, CFUN, CPLX, REAL, EPS)                 \
static int NAME(int is_host, double *rel_out, double *tol_out)             \
{                                                                         \
    STRUC id; memset(&id, 0, sizeof id);                                  \
    MUMPS_INT irn[N*N], jcn[N*N];                                         \
    CPLX a[N*N], rhs[N], xt[N], A[N][N];                                  \
    int i, j, k = 0;                                                      \
    xt[0].r=(REAL)1.0;  xt[0].i=(REAL)1.0;                                \
    xt[1].r=(REAL)-2.0; xt[1].i=(REAL)0.5;                                \
    xt[2].r=(REAL)3.0;  xt[2].i=(REAL)-1.0;                               \
    xt[3].r=(REAL)-4.0; xt[3].i=(REAL)2.0;                                \
    for (i=0;i<N;i++) for (j=0;j<N;j++){ A[i][j].r=(REAL)Ar[i][j]; A[i][j].i=(REAL)Ai[i][j]; } \
    for (j=0;j<N;j++)                                                     \
        for (i=0;i<N;i++,k++){ irn[k]=i+1; jcn[k]=j+1; a[k]=A[i][j]; }     \
    for (i=0;i<N;i++){                                                    \
        REAL sr=(REAL)0, si=(REAL)0;                                      \
        for (j=0;j<N;j++){                                               \
            sr += A[i][j].r*xt[j].r - A[i][j].i*xt[j].i;                  \
            si += A[i][j].r*xt[j].i + A[i][j].i*xt[j].r;                  \
        }                                                                \
        rhs[i].r=sr; rhs[i].i=si;                                         \
    }                                                                    \
    id.par = 1; id.sym = 0;                                               \
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);                       \
    id.job = -1; CFUN(&id);                                               \
    if (id.infog[0] < 0) return id.infog[0];                             \
    id.icntl[0]=-1; id.icntl[1]=-1; id.icntl[2]=-1; id.icntl[3]=0;       \
    if (is_host) { id.n=N; id.nnz=(MUMPS_INT8)(N*N); id.irn=irn; id.jcn=jcn; id.a=a; id.rhs=rhs; } \
    id.job = 6; CFUN(&id);                                                \
    if (id.infog[0] < 0) return id.infog[0];                             \
    if (is_host) {                                                        \
        double mr=0.0, den=0.0;                                           \
        for (i=0;i<N;i++){ double ax=(double)xt[i].r*xt[i].r+(double)xt[i].i*xt[i].i; if(ax>den) den=ax; } \
        for (i=0;i<N;i++){ double dr=(double)rhs[i].r-xt[i].r, di=(double)rhs[i].i-xt[i].i; \
                           double d=dr*dr+di*di; if(d>mr) mr=d; }         \
        double rel=(den>0.0)?mr/den:mr;                                   \
        *rel_out=sqrt(rel); *tol_out=16.0*(double)(N*N*N)*(EPS);          \
    }                                                                    \
    id.job = -2; CFUN(&id);                                               \
    return 0;                                                             \
}

GEN_REAL_SOLVE(solve_s, SMUMPS_STRUC_C, smumps_c, float,  (double)FLT_EPSILON)
GEN_REAL_SOLVE(solve_d, DMUMPS_STRUC_C, dmumps_c, double, DBL_EPSILON)
GEN_CPLX_SOLVE(solve_c, CMUMPS_STRUC_C, cmumps_c, mumps_complex,        float,  (double)FLT_EPSILON)
GEN_CPLX_SOLVE(solve_z, ZMUMPS_STRUC_C, zmumps_c, mumps_double_complex, double, DBL_EPSILON)

int main(int argc, char **argv)
{
    int myid, is_host;
    double rel = 0.0, tol = 0.0;
    struct { const char *label; int (*fn)(int, double *, double *); } cases[] = {
        {"single-real",    solve_s},
        {"single-complex", solve_c},
        {"double-real",    solve_d},
        {"double-complex", solve_z},
    };
    size_t nc = sizeof cases / sizeof cases[0], ci;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    is_host = (myid == 0);

    if (is_host) report_init("test_genuine_scdz_solve", "genuine");

    for (ci = 0; ci < nc; ci++) {
        int code = cases[ci].fn(is_host, &rel, &tol);
        if (code < 0) {
            fprintf(stderr, "genuine %s solve failed, infog=%d\n", cases[ci].label, code);
            if (is_host && gJson) report_finalize();
            MPI_Finalize();
            return 1;
        }
        if (is_host) report_case(cases[ci].label, rel, tol);
    }

    if (is_host) report_finalize();
    MPI_Finalize();
    return is_host ? (gAnyFail ? 1 : 0) : 0;
}
