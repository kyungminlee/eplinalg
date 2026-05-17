/* GENERATED-BY-gen_perf_harnesses — do not edit by hand; regenerate via
 *   python3 scripts/gen_perf_harnesses.py
 *
 * Kernel-isolated C perf harness for wtbmv (overlay vs migrated).
 * Built per-executable with -ffunction-sections / --gc-sections.
 */
#include "../perf_common.h"

#ifdef __cplusplus
#define BLAS_EXTERN extern "C"
#else
#define BLAS_EXTERN extern
#endif
typedef struct { double v[2]; } MFR;     /* float64x2 layout (POD) */
typedef struct { MFR r; MFR i; } MFC;    /* complex64x2 layout (POD) */
static inline MFR MFR_FROM(double d) { MFR x; x.v[0] = d; x.v[1] = 0.0; return x; }
static inline MFC MFC_FROM(double re, double im) {
    MFC z;
    z.r.v[0] = re; z.r.v[1] = 0.0;
    z.i.v[0] = im; z.i.v[1] = 0.0;
    return z;
}
static inline MFR Tr_from_d(double d) { return MFR_FROM(d); }
static inline MFC Tc_from_d(double d) { return MFC_FROM(d, 0.0); }


BLAS_EXTERN void wtbmv_(const char *, const char *, const char *, const int *, const int *,
    const MFC *, const int *, MFC *, const int *, size_t, size_t, size_t);
BLAS_EXTERN void wtbmv_migrated_(const char *, const char *, const char *, const int *, const int *,
    const MFC *, const int *, MFC *, const int *, size_t, size_t, size_t);

static void run_one(char uplo, char trans, char diag, int N, int K, int iters, int warmup) {
    int one = 1;
    int LDA = K + 1;
    MFC *A  = (MFC *)perf_aligned_alloc(64, (size_t)LDA * (size_t)N * sizeof(MFC));
    MFC *X  = (MFC *)perf_aligned_alloc(64, (size_t)N * sizeof(MFC));
    MFC *Xi = (MFC *)perf_aligned_alloc(64, (size_t)N * sizeof(MFC));
    for (size_t i = 0; i < (size_t)LDA*N; ++i) { int s = 2; A[i] = MFC_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    /* diagonal at known row of band — large to stabilize tbsv */
    int diag_row = (uplo == 'U') ? K : 0;
    for (int j = 0; j < N; ++j) A[(size_t)j * LDA + diag_row] = Tc_from_d((double)(K + 4));
    for (int i = 0; i < N; ++i) { int s = 3; Xi[i] = MFC_FROM(perf_fill_double(i, s), perf_fill_double(i, s + 131)); }
    memcpy(X, Xi, (size_t)N * sizeof(MFC));
    for (int r = 0; r < warmup; ++r) {
        wtbmv_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(MFC));
        wtbmv_migrated_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(MFC));
    }
    double t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        wtbmv_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(MFC));
    }
    double t1 = perf_now_s();
    double t_ov = (t1 - t0) / (iters ? iters : 1);
    t0 = perf_now_s();
    for (int it = 0; it < iters; ++it) {
        wtbmv_migrated_(&uplo, &trans, &diag, &N, &K, A, &LDA, X, &one, 1, 1, 1);
        memcpy(X, Xi, (size_t)N * sizeof(MFC));
    }
    t1 = perf_now_s();
    double t_mg = (t1 - t0) / (iters ? iters : 1);
    double flops = 4.0 * (double)(2*K+1) * (double)N;
    char key[4] = {uplo, trans, diag, 0};
    perf_emit("wtbmv", key, N, iters, flops, t_ov, t_mg);
    perf_emit_json("wtbmv", key, N, iters, flops, t_ov, t_mg);
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
    const char transes[] = { 'N','T','C' };
    for (size_t u = 0; u < 2; ++u) for (size_t t = 0; t < sizeof(transes); ++t) {
        char uplo = (u == 0) ? 'U' : 'L';
        char trans = transes[t];
        char diag = 'N';
        for (int i = 0; i < n; ++i) run_one(uplo, trans, diag, sizes[i], 16, iters, warmup);
    }
    return 0;
}
