/* perf_egemv — C perf harness for kind10 egemv overlay vs migrated.
 *
 * Why a C harness alongside the Fortran bench:
 *   The fypp-generated Fortran bench resets Y from Y_init each iter via
 *   a per-iter `realloc()` + scalar x87 copy loop. That introduces a
 *   measurement asymmetry that systematically penalizes the overlay by
 *   ~10% — not a kernel issue. This harness calls each routine in a
 *   tight isolated loop (no Y reset between calls) so the reported
 *   GF/s reflects steady-state kernel throughput.
 *
 * Env:
 *   BLAS_PERF_SIZES   comma-separated list, default "128,256,512,1024,2048"
 *   BLAS_PERF_ITERS   timed iters per (trans,size), default 30
 *   BLAS_PERF_WARMUP  warmup calls before timing, default 5
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef long double T;

extern void egemv_(const char *trans, const int *m, const int *n,
                   const T *alpha, const T *a, const int *lda,
                   const T *x, const int *incx,
                   const T *beta, T *y, const int *incy, size_t trans_len);

extern void egemv_migrated_(const char *trans, const int *m, const int *n,
                            const T *alpha, const T *a, const int *lda,
                            const T *x, const int *incx,
                            const T *beta, T *y, const int *incy, size_t trans_len);

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}

static int parse_sizes(int *out, int max) {
    const char *s = getenv("BLAS_PERF_SIZES");
    if (!s || !*s) {
        const int defaults[] = {128, 256, 512, 1024, 2048};
        int n = (int)(sizeof(defaults)/sizeof(defaults[0]));
        if (n > max) n = max;
        for (int i = 0; i < n; ++i) out[i] = defaults[i];
        return n;
    }
    int n = 0;
    const char *p = s;
    while (*p && n < max) {
        char *end;
        long v = strtol(p, &end, 10);
        if (end == p) break;
        if (v > 0) out[n++] = (int)v;
        p = end;
        while (*p == ',' || *p == ' ') ++p;
    }
    return n;
}

static void run_one(char trans, int N, int iters, int warmup) {
    int M = N;
    T alpha = 0.7L, beta = 0.3L;
    int one = 1;

    T *A  = aligned_alloc(64, (size_t)M * N * sizeof(T));
    T *X  = aligned_alloc(64, (size_t)N * sizeof(T));
    T *Y  = aligned_alloc(64, (size_t)M * sizeof(T));
    if (!A || !X || !Y) { fprintf(stderr, "alloc failed\n"); exit(2); }

    for (size_t i = 0; i < (size_t)M * N; ++i) A[i] = (T)((i * 31 + 7) % 17) / 13.0L;
    for (int i = 0; i < N; ++i) X[i] = (T)((i * 5) % 11) / 7.0L;
    for (int i = 0; i < M; ++i) Y[i] = (T)((i * 3) % 13) / 11.0L;

    for (int r = 0; r < warmup; ++r) {
        egemv_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
        egemv_migrated_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
    }

    /* Block timing: one clock pair per `iters` calls. Avoids per-call
     * clock overhead and amortizes any first-call cache-warm effects.
     * Reset Y to a known-bounded vector once before each block so the
     * compounding Y = beta*Y + alpha*A*x updates don't drift into
     * x87-slow ranges (denormals, very large magnitudes). */
    T *Yi = aligned_alloc(64, (size_t)M * sizeof(T));
    for (int i = 0; i < M; ++i) Yi[i] = (T)((i * 3) % 13) / 11.0L;

    double t0, t1;
    memcpy(Y, Yi, (size_t)M * sizeof(T));
    t0 = now_s();
    for (int it = 0; it < iters; ++it)
        egemv_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
    t1 = now_s();
    double t_ov = (t1 - t0) / iters;

    memcpy(Y, Yi, (size_t)M * sizeof(T));
    t0 = now_s();
    for (int it = 0; it < iters; ++it)
        egemv_migrated_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
    t1 = now_s();
    double t_mg = (t1 - t0) / iters;

    free(Yi);
    double flops = 2.0 * (double)M * (double)N;
    double g_ov = (t_ov > 0) ? flops / t_ov / 1e9 : 0;
    double g_mg = (t_mg > 0) ? flops / t_mg / 1e9 : 0;
    double ratio = (t_ov > 0) ? t_mg / t_ov : 0;

    printf("%c   %7d   %12.4f   %12.4f   %8.3fx   sink=%.3Lg\n",
           trans, N, g_ov, g_mg, ratio, Y[0]);

    free(A); free(X); free(Y);
}

int main(void) {
    int iters  = env_int("BLAS_PERF_ITERS",  30);
    int warmup = env_int("BLAS_PERF_WARMUP", 5);
    int sizes[32];
    int nsizes = parse_sizes(sizes, 32);

    printf("trans    size      overlay GF/s   migrated GF/s   mig/ov\n");
    const char transes[] = {'N', 'T'};
    for (size_t t = 0; t < sizeof(transes); ++t) {
        for (int i = 0; i < nsizes; ++i) {
            run_one(transes[t], sizes[i], iters, warmup);
        }
    }
    return 0;
}
