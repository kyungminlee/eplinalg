/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for qger (overlay vs migrated).
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


BLAS_EXTERN void qger_(const int *, const int *, const Q16 *,
    const Q16 *, const int *, const Q16 *, const int *, Q16 *, const int *);
BLAS_EXTERN void qger_migrated_(const int *, const int *, const Q16 *,
    const Q16 *, const int *, const Q16 *, const int *, Q16 *, const int *);

static void run_one(int M, int N, int iters, int warmup) {
    int one = 1;
    Q16 alpha = Q16_FROM(0.7);
    Q16 *A  = (Q16 *)perf_aligned_alloc(64, (size_t)M * (size_t)N * sizeof(Q16));
    Q16 *Ai = (Q16 *)perf_aligned_alloc(64, (size_t)M * (size_t)N * sizeof(Q16));
    Q16 *X  = (Q16 *)perf_aligned_alloc(64, (size_t)M * sizeof(Q16));
    Q16 *Y  = (Q16 *)perf_aligned_alloc(64, (size_t)N * sizeof(Q16));
    for (size_t i = 0; i < (size_t)M*N; ++i) { int s = 2; Ai[i] = Q16_FROM(perf_fill_double(i, s)); }
    for (int i = 0; i < M; ++i)           { int s = 3; X[i] = Q16_FROM(perf_fill_double(i, s)); }
    for (int i = 0; i < N; ++i)           { int s = 4; Y[i] = Q16_FROM(perf_fill_double(i, s)); }
    memcpy(A, Ai, (size_t)M * (size_t)N * sizeof(Q16));

    for (int r = 0; r < warmup; ++r) {
        qger_(&M, &N, &alpha, X, &one, Y, &one, A, &M);
        memcpy(A, Ai, (size_t)M * (size_t)N * sizeof(Q16));
        qger_migrated_(&M, &N, &alpha, X, &one, Y, &one, A, &M);
        memcpy(A, Ai, (size_t)M * (size_t)N * sizeof(Q16));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        qger_(&M, &N, &alpha, X, &one, Y, &one, A, &M);
        memcpy(A, Ai, (size_t)M * (size_t)N * sizeof(Q16));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);

    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        qger_migrated_(&M, &N, &alpha, X, &one, Y, &one, A, &M);
        memcpy(A, Ai, (size_t)M * (size_t)N * sizeof(Q16));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);

    double flops = 2.0 * (double)M * (double)N;
    perf_emit("qger", "-", N, iters, flops, t_ov, t_mg);
    perf_emit_json("qger", "-", N, iters, flops, t_ov, t_mg);
    free(A); free(Ai); free(X); free(Y);
}

static const int default_sizes[] = {128, 256, 512, 1024};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  200);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 20);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    for (int i = 0; i < n; ++i) run_one(sizes[i], sizes[i], iters, warmup);
    return 0;
}
