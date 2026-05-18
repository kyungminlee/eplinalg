/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for ygemv (overlay vs migrated).
 * Built per-executable with -ffunction-sections / --gc-sections.
 */
#include "../perf_common.h"

#include <complex.h>
#ifdef __cplusplus
#define BLAS_EXTERN extern "C"
#else
#define BLAS_EXTERN extern
#endif
typedef long double R10;
typedef _Complex long double C10;
#define R10_FROM(d) ((R10)(d))
#define C10_FROM(re, im) ((R10)(re) + 1.0iL * (R10)(im))
static inline R10 Tr_from_d(double d) { return (R10)d; }
static inline C10 Tc_from_d(double d) { return (C10)d; }


BLAS_EXTERN void ygemv_(const char *, const int *, const int *, const C10 *,
    const C10 *, const int *, const C10 *, const int *,
    const C10 *, C10 *, const int *, size_t);
BLAS_EXTERN void ygemv_migrated_(const char *, const int *, const int *, const C10 *,
    const C10 *, const int *, const C10 *, const int *,
    const C10 *, C10 *, const int *, size_t);

static void run_one(char trans, int M, int N, int iters, int warmup) {
    int one = 1;
    C10 alpha = C10_FROM(0.7, 0.0), beta = C10_FROM(0.3, 0.0);
    C10 *A  = (C10 *)perf_aligned_alloc(64, (size_t)M * (size_t)N * sizeof(C10));
    C10 *X  = (C10 *)perf_aligned_alloc(64, (size_t)N * sizeof(C10));
    C10 *Y  = (C10 *)perf_aligned_alloc(64, (size_t)M * sizeof(C10));
    C10 *Yi = (C10 *)perf_aligned_alloc(64, (size_t)M * sizeof(C10));
    for (size_t i = 0; i < (size_t)M*N; ++i) { int s = 2; A[i] = C10_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    for (int i = 0; i < N; ++i)           { int s = 3; X[i] = C10_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    for (int i = 0; i < M; ++i)           { int s = 4; Yi[i] = C10_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }

    memcpy(Y, Yi, (size_t)M * sizeof(C10));
    for (int r = 0; r < warmup; ++r) {
        ygemv_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
        memcpy(Y, Yi, (size_t)M * sizeof(C10));
        ygemv_migrated_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
        memcpy(Y, Yi, (size_t)M * sizeof(C10));
    }

    memcpy(Y, Yi, (size_t)M * sizeof(C10));
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it)
        ygemv_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);

    memcpy(Y, Yi, (size_t)M * sizeof(C10));
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it)
        ygemv_migrated_(&trans, &M, &N, &alpha, A, &M, X, &one, &beta, Y, &one, 1);
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);

    double flops = 8.0 * (double)M * (double)N;
    char key[2] = {trans, 0};
    perf_emit("ygemv", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("ygemv", key, N, iters, flops, t_ov, t_mg);
    free(A); free(X); free(Y); free(Yi);
}

static const int default_sizes[] = {128, 256, 512, 1024, 2048};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  200);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 20);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    const char transes[] = { 'N','T','C' };
    for (size_t t = 0; t < sizeof(transes); ++t)
        for (int i = 0; i < n; ++i)
            run_one(transes[t], sizes[i], sizes[i], iters, warmup);
    return 0;
}
