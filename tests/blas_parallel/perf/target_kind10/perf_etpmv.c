/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for etpmv (overlay vs migrated).
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


BLAS_EXTERN void etpmv_(const char *, const char *, const char *, const int *,
    const R10 *, R10 *, const int *, size_t, size_t, size_t);
BLAS_EXTERN void etpmv_migrated_(const char *, const char *, const char *, const int *,
    const R10 *, R10 *, const int *, size_t, size_t, size_t);

static void run_one(char uplo, char trans, char diag, int N, int iters, int warmup) {
    int one = 1;
    size_t AP_LEN = (size_t)N * (size_t)(N + 1) / 2;
    R10 *AP = (R10 *)perf_aligned_alloc(64, AP_LEN * sizeof(R10));
    R10 *X  = (R10 *)perf_aligned_alloc(64, (size_t)N * sizeof(R10));
    R10 *Xi = (R10 *)perf_aligned_alloc(64, (size_t)N * sizeof(R10));
    for (size_t i = 0; i < AP_LEN; ++i) { int s = 2; AP[i] = R10_FROM(perf_fill_double(i, s)); }
    /* Force diagonal to ~N for stability of tpsv */
    if (uplo == 'U') {
        size_t off = 0;
        for (int j = 0; j < N; ++j) { AP[off + j] = Tr_from_d((double)(N + 4)); off += (size_t)(j + 1); }
    } else {
        size_t off = 0;
        for (int j = 0; j < N; ++j) { AP[off] = Tr_from_d((double)(N + 4)); off += (size_t)(N - j); }
    }
    for (int i = 0; i < N; ++i) { int s = 3; Xi[i] = R10_FROM(perf_fill_double(i, s)); }
    memcpy(X, Xi, (size_t)N * sizeof(R10));
    for (int r = 0; r < warmup; ++r) {
        etpmv_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
        etpmv_migrated_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        etpmv_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        etpmv_migrated_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 1.0 * (double)N * (double)N;
    char key[4] = {uplo, trans, diag, 0};
    perf_emit("etpmv", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("etpmv", key, N, iters, flops, t_ov, t_mg);
    free(AP); free(X); free(Xi);
}

static const int default_sizes[] = {128, 256, 512, 1024};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  20);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 3);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    const char transes[] = { 'N','T' };
    for (size_t u = 0; u < 2; ++u) for (size_t t = 0; t < sizeof(transes); ++t) {
        char uplo = (u == 0) ? 'U' : 'L';
        char trans = transes[t];
        char diag = 'N';
        for (int i = 0; i < n; ++i) run_one(uplo, trans, diag, sizes[i], iters, warmup);
    }
    return 0;
}
