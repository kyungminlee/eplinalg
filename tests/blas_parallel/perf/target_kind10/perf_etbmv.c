/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for etbmv (overlay vs migrated).
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


BLAS_EXTERN void etbmv_(const char *, const char *, const char *, const int *, const int *,
    const R10 *, const int *, R10 *, const int *, size_t, size_t, size_t);
BLAS_EXTERN void etbmv_migrated_(const char *, const char *, const char *, const int *, const int *,
    const R10 *, const int *, R10 *, const int *, size_t, size_t, size_t);

static void run_one(char uplo, char trans, char diag, int N, int K, int iters, int warmup) {
    int one = 1;
    int LDA = K + 1;
    R10 *A  = (R10 *)perf_aligned_alloc(64, (size_t)LDA * (size_t)N * sizeof(R10));
    R10 *X  = (R10 *)perf_aligned_alloc(64, (size_t)N * sizeof(R10));
    R10 *Xi = (R10 *)perf_aligned_alloc(64, (size_t)N * sizeof(R10));
    for (size_t i = 0; i < (size_t)LDA*N; ++i) { int s = 2; A[i] = R10_FROM(perf_fill_double(i, s)); }
    /* diagonal at known row of band — large to stabilize tbsv */
    int diag_row = (uplo == 'U') ? K : 0;
    for (int j = 0; j < N; ++j) A[(size_t)j * LDA + diag_row] = Tr_from_d((double)(K + 4));
    for (int i = 0; i < N; ++i) { int s = 3; Xi[i] = R10_FROM(perf_fill_double(i, s)); }
    memcpy(X, Xi, (size_t)N * sizeof(R10));
    for (int r = 0; r < warmup; ++r) {
        etbmv_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
        etbmv_migrated_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        etbmv_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        etbmv_migrated_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(R10));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 1.0 * (double)(2*K+1) * (double)N;
    char key[4] = {uplo, trans, diag, 0};
    perf_emit("etbmv", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("etbmv", key, N, iters, flops, t_ov, t_mg);
    free(A); free(X); free(Xi);
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
        for (int i = 0; i < n; ++i) run_one(uplo, trans, diag, sizes[i], 16, iters, warmup);
    }
    return 0;
}
