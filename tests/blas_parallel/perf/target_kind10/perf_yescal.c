/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for yescal (overlay vs migrated).
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


BLAS_EXTERN void yescal_(const int *, const R10 *, C10 *, const int *);
BLAS_EXTERN void yescal_migrated_(const int *, const R10 *, C10 *, const int *);

static void run_yescal(int N, int iters, int warmup) {
    int one = 1;
    R10 alpha = R10_FROM(0.7);
    C10 *X = (C10 *)perf_aligned_alloc(64, (size_t)N * sizeof(C10));
    C10 *Xi = (C10 *)perf_aligned_alloc(64, (size_t)N * sizeof(C10));
    for (int i = 0; i < N; ++i) { int s = 0; Xi[i] = C10_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    memcpy(X, Xi, (size_t)N * sizeof(C10));
    for (int r = 0; r < warmup; ++r) {
        yescal_(&N, &alpha, X, &one);
        memcpy(X, Xi, (size_t)N * sizeof(C10));
        yescal_migrated_(&N, &alpha, X, &one);
        memcpy(X, Xi, (size_t)N * sizeof(C10));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        yescal_(&N, &alpha, X, &one);
        memcpy(X, Xi, (size_t)N * sizeof(C10));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);

    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        yescal_migrated_(&N, &alpha, X, &one);
        memcpy(X, Xi, (size_t)N * sizeof(C10));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 6.0 * (double)N;
    perf_emit("yescal", "-", N, iters, flops, t_ov, t_mg);
    perf_emit_json("yescal", "-", N, iters, flops, t_ov, t_mg);
    free(X); free(Xi);
}

static const int default_sizes[] = {4096, 16384, 65536};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  20);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 3);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    for (int i = 0; i < n; ++i) run_yescal(sizes[i], iters, warmup);
    return 0;
}
