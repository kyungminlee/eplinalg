/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for xtrmm (overlay vs migrated).
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


BLAS_EXTERN void xtrmm_(const char *, const char *, const char *, const char *,
    const int *, const int *, const X16 *,
    const X16 *, const int *, X16 *, const int *,
    size_t, size_t, size_t, size_t);
BLAS_EXTERN void xtrmm_migrated_(const char *, const char *, const char *, const char *,
    const int *, const int *, const X16 *,
    const X16 *, const int *, X16 *, const int *,
    size_t, size_t, size_t, size_t);

static void run_one(char side, char uplo, char trans, char diag,
                    int M, int N, int iters, int warmup) {
    X16 alpha = X16_FROM(0.7, 0.0);
    int Asz = (side == 'L') ? M : N;
    int lda = Asz, ldb = M;
    X16 *A  = (X16 *)perf_aligned_alloc(64, (size_t)Asz * (size_t)Asz * sizeof(X16));
    X16 *B  = (X16 *)perf_aligned_alloc(64, (size_t)M * (size_t)N * sizeof(X16));
    X16 *Bi = (X16 *)perf_aligned_alloc(64, (size_t)M * (size_t)N * sizeof(X16));
    for (size_t i = 0; i < (size_t)Asz*Asz; ++i) { int s = 2; A[i] = X16_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    /* diagonal dominance for trsm */
    for (int i = 0; i < Asz; ++i) A[(size_t)i * lda + i] = Tc_from_d((double)(Asz + 4));
    for (size_t i = 0; i < (size_t)M*N; ++i) { int s = 4; Bi[i] = X16_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    memcpy(B, Bi, (size_t)M * (size_t)N * sizeof(X16));
    for (int r = 0; r < warmup; ++r) {
        xtrmm_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
        memcpy(B, Bi, (size_t)M * (size_t)N * sizeof(X16));
        xtrmm_migrated_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
        memcpy(B, Bi, (size_t)M * (size_t)N * sizeof(X16));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        xtrmm_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
        memcpy(B, Bi, (size_t)M * (size_t)N * sizeof(X16));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        xtrmm_migrated_(&side, &uplo, &trans, &diag, &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
        memcpy(B, Bi, (size_t)M * (size_t)N * sizeof(X16));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 4.0 * (double)M * (double)N * (double)M;
    char key[5] = {side, uplo, trans, diag, 0};
    perf_emit("xtrmm", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("xtrmm", key, N, iters, flops, t_ov, t_mg);
    free(A); free(B); free(Bi);
}

static const int default_sizes[] = {64, 128, 256, 512};
int main(void) {
    int iters  = perf_env_int("BLAS_PERF_ITERS",  10);
    int warmup = perf_env_int("BLAS_PERF_WARMUP", 2);
    int sizes[32];
    int n = perf_parse_sizes(default_sizes,
        (int)(sizeof(default_sizes)/sizeof(default_sizes[0])), sizes, 32);
    perf_print_header();
    /* Sample a small set of (side, uplo, trans, diag) — not full 16. */
    const char sides[] = {'L', 'R'};
    const char uplos[] = {'U', 'L'};
    const char transes[] = { 'N', 'T', 'C' };
    for (size_t s = 0; s < 2; ++s) for (size_t u = 0; u < 2; ++u)
      for (size_t t = 0; t < sizeof(transes); ++t)
        for (int i = 0; i < n; ++i)
            run_one(sides[s], uplos[u], transes[t], 'N', sizes[i], sizes[i], iters, warmup);
    return 0;
}
