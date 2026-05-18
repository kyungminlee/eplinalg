/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for xtpmv (overlay vs migrated).
 * Built per-executable with -ffunction-sections / --gc-sections.
 */
#include "../perf_common.h"

#include <quadmath.h>
#ifdef __cplusplus
#define BLAS_EXTERN extern "C"
#else
#define BLAS_EXTERN extern
#endif
typedef __float128 Q16;
typedef _Complex float __attribute__((mode(TC))) X16;
#define Q16_FROM(d) ((Q16)(double)(d))
#define X16_FROM(re, im) ((X16)((Q16)(double)(re) + 1.0i * (Q16)(double)(im)))
static inline Q16 Tr_from_d(double d) { return (Q16)d; }
static inline X16 Tc_from_d(double d) { return (X16)((Q16)d); }


BLAS_EXTERN void xtpmv_(const char *, const char *, const char *, const int *,
    const X16 *, X16 *, const int *, size_t, size_t, size_t);
BLAS_EXTERN void xtpmv_migrated_(const char *, const char *, const char *, const int *,
    const X16 *, X16 *, const int *, size_t, size_t, size_t);

static void run_one(char uplo, char trans, char diag, int N, int iters, int warmup) {
    int one = 1;
    size_t AP_LEN = (size_t)N * (size_t)(N + 1) / 2;
    X16 *AP = (X16 *)perf_aligned_alloc(64, AP_LEN * sizeof(X16));
    X16 *X  = (X16 *)perf_aligned_alloc(64, (size_t)N * sizeof(X16));
    X16 *Xi = (X16 *)perf_aligned_alloc(64, (size_t)N * sizeof(X16));
    for (size_t i = 0; i < AP_LEN; ++i) { int s = 2; AP[i] = X16_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    /* Force diagonal to ~N for stability of tpsv */
    if (uplo == 'U') {
        size_t off = 0;
        for (int j = 0; j < N; ++j) { AP[off + j] = Tc_from_d((double)(N + 4)); off += (size_t)(j + 1); }
    } else {
        size_t off = 0;
        for (int j = 0; j < N; ++j) { AP[off] = Tc_from_d((double)(N + 4)); off += (size_t)(N - j); }
    }
    for (int i = 0; i < N; ++i) { int s = 3; Xi[i] = X16_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    memcpy(X, Xi, (size_t)N * sizeof(X16));
    for (int r = 0; r < warmup; ++r) {
        xtpmv_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(X16));
        xtpmv_migrated_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(X16));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        xtpmv_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(X16));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        xtpmv_migrated_(&uplo, &trans, &diag, &N, AP, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(X16));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 4.0 * (double)N * (double)N;
    char key[4] = {uplo, trans, diag, 0};
    perf_emit("xtpmv", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("xtpmv", key, N, iters, flops, t_ov, t_mg);
    free(AP); free(X); free(Xi);
}

static const int default_sizes[] = {128, 256, 512, 1024};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  200);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 20);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    const char transes[] = { 'N','T','C' };
    for (size_t u = 0; u < 2; ++u) for (size_t t = 0; t < sizeof(transes); ++t) {
        char uplo = (u == 0) ? 'U' : 'L';
        char trans = transes[t];
        char diag = 'N';
        for (int i = 0; i < n; ++i) run_one(uplo, trans, diag, sizes[i], iters, warmup);
    }
    return 0;
}
