/*
 * egemm — kind10 (REAL(KIND=10), x86-64 80-bit long double) GEMM overlay.
 *
 * Step 3: cache-blocked + packed implementation for the no-transpose
 * NN case (the LAPACK-dominant code path). The other 8 trans combos
 * fall through to the naive triple-loop reference. Step 4 (next)
 * adds OpenMP across the NC loop.
 *
 * Why only cache blocking (no register tiling) on kind10:
 *   - x86-64 `long double` lives in the x87 stack — 8 registers, all
 *     of which see frequent spills. A 4x4 register tile that works
 *     beautifully for `double` (16 xmm registers, FMA, AVX2) buys
 *     little here. The dominant win for kind10 is eliminating L1/L2
 *     stride misses via packing.
 *   - We still keep a tight inner loop where the compiler can hoist
 *     the alpha*B(p,j) factor; that gives us roughly the same shape
 *     as a 1xNR register tile.
 *
 * Layout:
 *   A panel: packed column-major, contiguous within each KC-row band,
 *     contiguous across rows of each MC block. Read sequentially in
 *     the inner kernel.
 *   B block: packed column-major (KC rows × NC cols), so for each
 *     fixed j the K-axis traversal of B is contiguous.
 *
 * Block sizes (tunable via env at step 3.5 autotune):
 *   MC * KC fits in L2  (kind10: 16 bytes/elt → 64*128 = 128 KiB).
 *   KC * NC fits in L3  (rarely the constraint at the sizes we target).
 *
 * Fortran ABI:
 *   - subroutine name lowercased + trailing underscore: `egemm_`
 *   - scalars passed by pointer
 *   - character args followed by hidden trailing `size_t` lengths
 *   - REAL(KIND=10) ↔ `long double` (x86-64 80-bit extended)
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef long double T;

/* Block sizes. Env-overridable for autotune (step 3.5). */
static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}

static int g_mc = 0, g_kc = 0, g_nc = 0;
static void init_blocks(void) {
    if (g_mc) return;
    g_mc = env_int("EBLAS_MC",  64);
    g_kc = env_int("EBLAS_KC", 128);
    g_nc = env_int("EBLAS_NC", 256);
}

static int trans_code(const char *p, size_t len) {
    (void)len;
    char c = (char)toupper((unsigned char)*p);
    return c;
}

/* ── Naive reference for non-NN paths and tiny problems ────────── */

static void naive_nn(int m, int n, int k, T alpha,
                     const T *a, int lda, const T *b, int ldb,
                     T beta, T *c, int ldc)
{
    int i, j, p;
    for (j = 0; j < n; ++j) {
        T *cj = &c[(size_t)j * ldc];
        if (beta == 0.0L)      for (i = 0; i < m; ++i) cj[i]  = 0.0L;
        else if (beta != 1.0L) for (i = 0; i < m; ++i) cj[i] *= beta;
        const T *bj = &b[(size_t)j * ldb];
        for (p = 0; p < k; ++p) {
            const T t = alpha * bj[p];
            const T *ap = &a[(size_t)p * lda];
            for (i = 0; i < m; ++i) cj[i] += t * ap[i];
        }
    }
}

static void naive_other(int ta, int tb, int m, int n, int k, T alpha,
                        const T *a, int lda, const T *b, int ldb,
                        T beta, T *c, int ldc)
{
    int i, j, p;
    for (j = 0; j < n; ++j) {
        T *cj = &c[(size_t)j * ldc];
        if (beta == 0.0L)      for (i = 0; i < m; ++i) cj[i]  = 0.0L;
        else if (beta != 1.0L) for (i = 0; i < m; ++i) cj[i] *= beta;
    }
    if (ta == 'N' && tb != 'N') {
        for (j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            for (p = 0; p < k; ++p) {
                const T t = alpha * b[(size_t)p * ldb + j];
                const T *ap = &a[(size_t)p * lda];
                for (i = 0; i < m; ++i) cj[i] += t * ap[i];
            }
        }
    } else if (ta != 'N' && tb == 'N') {
        for (j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            const T *bj = &b[(size_t)j * ldb];
            for (i = 0; i < m; ++i) {
                const T *ai = &a[(size_t)i * lda];
                T acc = 0.0L;
                for (p = 0; p < k; ++p) acc += ai[p] * bj[p];
                cj[i] += alpha * acc;
            }
        }
    } else {
        for (j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            for (i = 0; i < m; ++i) {
                const T *ai = &a[(size_t)i * lda];
                T acc = 0.0L;
                for (p = 0; p < k; ++p) acc += ai[p] * b[(size_t)p * ldb + j];
                cj[i] += alpha * acc;
            }
        }
    }
}

/* ── Packed + blocked NN ──────────────────────────────────────── */

/*
 * Pack A(ic:ic+ib, pc:pc+pb) into Ap, column-major MR×KC stripes.
 * Inner kernel reads Ap[p*ib + i], i.e. column-major (ib, pb). We
 * keep it that simple — no MR-row stripes — because the inner loop
 * over i is already cache-line-friendly after packing.
 */
static void pack_A_nn(const T *restrict A, int lda,
                      int ic, int pc, int ib, int pb,
                      T *restrict Ap)
{
    int i, p;
    for (p = 0; p < pb; ++p) {
        const T *src = &A[(size_t)(pc + p) * lda + ic];
        T *dst = &Ap[(size_t)p * ib];
        for (i = 0; i < ib; ++i) dst[i] = src[i];
    }
}

/*
 * Pack B(pc:pc+pb, jc:jc+jb) into Bp, column-major (pb, jb).
 * Same layout as the source — copy is for L2 locality, not reshape.
 */
static void pack_B_nn(const T *restrict B, int ldb,
                      int pc, int jc, int pb, int jb,
                      T *restrict Bp)
{
    int p, j;
    for (j = 0; j < jb; ++j) {
        const T *src = &B[(size_t)(jc + j) * ldb + pc];
        T *dst = &Bp[(size_t)j * pb];
        for (p = 0; p < pb; ++p) dst[p] = src[p];
    }
}

/*
 * C(ic:ic+ib, jc:jc+jb) += alpha * Ap[ib, pb] * Bp[pb, jb].
 * Ap is (ib × pb) column-major; Bp is (pb × jb) column-major.
 * The outer-j loop hoists alpha * Bp[p, j] into a scalar; the inner
 * i loop is one contiguous Ap stripe per p.
 */
static void inner_kernel_nn(int ib, int jb, int pb, T alpha,
                            const T *restrict Ap, const T *restrict Bp,
                            T *restrict C, int ldc)
{
    int i, j, p;
    for (j = 0; j < jb; ++j) {
        T *cj = &C[(size_t)j * ldc];
        const T *bj = &Bp[(size_t)j * pb];
        for (p = 0; p < pb; ++p) {
            const T t = alpha * bj[p];
            const T *ap = &Ap[(size_t)p * ib];
            for (i = 0; i < ib; ++i) cj[i] += t * ap[i];
        }
    }
}

static void blocked_nn(int M, int N, int K, T alpha,
                       const T *A, int lda, const T *B, int ldb,
                       T beta, T *C, int ldc)
{
    /* Apply beta to C up front. */
    int i, j;
    for (j = 0; j < N; ++j) {
        T *cj = &C[(size_t)j * ldc];
        if (beta == 0.0L)      for (i = 0; i < M; ++i) cj[i]  = 0.0L;
        else if (beta != 1.0L) for (i = 0; i < M; ++i) cj[i] *= beta;
    }
    if (alpha == 0.0L || K == 0) return;

    init_blocks();
    const int MC = g_mc, KC = g_kc, NC = g_nc;

    T *Ap = aligned_alloc(64, (size_t)MC * KC * sizeof(T));
    T *Bp = aligned_alloc(64, (size_t)KC * NC * sizeof(T));
    if (!Ap || !Bp) {
        /* aligned_alloc failed; fall back to naive. */
        free(Ap); free(Bp);
        naive_nn(M, N, K, alpha, A, lda, B, ldb, 1.0L, C, ldc);
        return;
    }

    int jc, pc, ic;
    for (jc = 0; jc < N; jc += NC) {
        const int jb = (N - jc < NC) ? (N - jc) : NC;
        for (pc = 0; pc < K; pc += KC) {
            const int pb = (K - pc < KC) ? (K - pc) : KC;
            pack_B_nn(B, ldb, pc, jc, pb, jb, Bp);
            for (ic = 0; ic < M; ic += MC) {
                const int ib = (M - ic < MC) ? (M - ic) : MC;
                pack_A_nn(A, lda, ic, pc, ib, pb, Ap);
                inner_kernel_nn(ib, jb, pb, alpha, Ap, Bp,
                                &C[(size_t)jc * ldc + ic], ldc);
            }
        }
    }

    free(Ap);
    free(Bp);
}

/* ── Entry point ──────────────────────────────────────────────── */

void egemm_(
    const char *transa, const char *transb,
    const int *m_, const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    size_t transa_len, size_t transb_len)
{
    const int m = *m_, n = *n_, k = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    if (m <= 0 || n <= 0) return;

    if (ta == 'N' && tb == 'N') {
        blocked_nn(m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    } else {
        naive_other(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
    }
}
