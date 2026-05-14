/*
 * egemm — kind10 (REAL(KIND=10), x86-64 80-bit long double) GEMM overlay.
 *
 * Goto-style cache-blocked + packed, OpenMP-parallel across the
 * outer NC loop. All 9 transpose combos go through the same blocked
 * path; the packers absorb TRANSA / TRANSB orientation.
 *
 * Inner-product (DDOT) micro-kernel: Ap is packed in ROW-MAJOR
 * layout (`Ap[i*pb + p] = op(A)[i,p]`), so the inner contraction
 * reads one contiguous Ap row and one contiguous Bp column as two
 * length-pb vectors and accumulates their dot product in a register.
 * C[i,j] is read+written once per (i,j) — half the C traffic of the
 * rank-1 outer-product variant. The accumulator can stay in a
 * register across the contraction.
 *
 * For real types (kind10) the conjugate-transpose 'C' is identical
 * to 'T'. We normalize both to 'T' at entry.
 *
 * No register tiling on kind10: x86-64 `long double` lives in the
 * 8-deep x87 stack — a 4×4 register tile spills immediately and
 * buys little. The dominant win comes from packing (kills L1/L2
 * stride misses on transposed inputs) and OpenMP across NC.
 *
 * Block sizes (env-overridable, autotune lands later):
 *   EBLAS_MC=64   panel rows
 *   EBLAS_KC=128  panel depth
 *   EBLAS_NC=256  column band per thread
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
#ifdef _OPENMP
#include <omp.h>
#endif

typedef long double T;

/* ── Block sizes ──────────────────────────────────────────────── */

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
    /* For real types, 'C' is identical to 'T'. */
    return (c == 'C') ? 'T' : c;
}

/* ── Packers ──────────────────────────────────────────────────── */

/*
 * Pack op(A)(ic:ic+ib, pc:pc+pb) into Ap in ROW-MAJOR layout
 * (`Ap[i*pb + p] = op(A)[i,p]`). Each row of Ap is contiguous along
 * the contraction axis p — what the inner-product kernel reads.
 *
 *   ta == 'N':  op(A)[i,p] = A[ic+i, pc+p].  Source stride in p is
 *               lda → gather across columns when writing one Ap row.
 *
 *   ta == 'T':  op(A)[i,p] = A[pc+p, ic+i].  Source stride in p is 1
 *               → straight memcpy per row of Ap.
 */
static void pack_A(const T *restrict A, int lda,
                   int ic, int pc, int ib, int pb,
                   int ta, T *restrict Ap)
{
    int i, p;
    if (ta == 'N') {
        for (i = 0; i < ib; ++i) {
            T *dst = &Ap[(size_t)i * pb];
            for (p = 0; p < pb; ++p)
                dst[p] = A[(size_t)(pc + p) * lda + (ic + i)];
        }
    } else {
        for (i = 0; i < ib; ++i) {
            const T *src = &A[(size_t)(ic + i) * lda + pc];
            T *dst = &Ap[(size_t)i * pb];
            for (p = 0; p < pb; ++p) dst[p] = src[p];
        }
    }
}

/*
 * Pack op(B)(pc:pc+pb, jc:jc+jb) into Bp, column-major (pb, jb).
 *
 *   tb == 'N':  op(B)[p, j] = B[pc+p, jc+j]
 *               source stride in p: 1
 *               source stride in j: ldb
 *
 *   tb == 'T':  op(B)[p, j] = B[jc+j, pc+p]
 *               source stride in p: ldb
 *               source stride in j: 1
 *
 * Pack output is always: Bp[j*pb + p] = op(B)[p, j].
 */
static void pack_B(const T *restrict B, int ldb,
                   int pc, int jc, int pb, int jb,
                   int tb, T *restrict Bp)
{
    int p, j;
    if (tb == 'N') {
        for (j = 0; j < jb; ++j) {
            const T *src = &B[(size_t)(jc + j) * ldb + pc];
            T *dst = &Bp[(size_t)j * pb];
            for (p = 0; p < pb; ++p) dst[p] = src[p];
        }
    } else {
        /* Walk row p of source B (stride ldb), scattering across
         * the j-axis of Bp. Outer p loop iterates ldb-stride rows. */
        for (p = 0; p < pb; ++p) {
            const T *src = &B[(size_t)(pc + p) * ldb + jc];
            for (j = 0; j < jb; ++j) Bp[(size_t)j * pb + p] = src[j];
        }
    }
}

/* ── Inner kernel ─────────────────────────────────────────────── */

/*
 * Inner-product micro-kernel.
 *   Ap row-major: Ap[i*pb + p] = op(A)[i, p]  (length-pb rows)
 *   Bp col-major: Bp[j*pb + p] = op(B)[p, j]  (length-pb cols)
 * For each (i, j) compute sum = Σ_p Ap[i,p] · Bp[p,j] in a scalar
 * accumulator, then C[i,j] += alpha · sum. One read+write of C per
 * (i,j), instead of one per (i,j,p).
 */
static void inner_kernel(int ib, int jb, int pb, T alpha,
                         const T *restrict Ap, const T *restrict Bp,
                         T *restrict C, int ldc)
{
    int i, j, p;
    for (j = 0; j < jb; ++j) {
        T *cj = &C[(size_t)j * ldc];
        const T *bj = &Bp[(size_t)j * pb];
        for (i = 0; i < ib; ++i) {
            const T *ai = &Ap[(size_t)i * pb];
            T sum = 0.0L;
            for (p = 0; p < pb; ++p) sum += ai[p] * bj[p];
            cj[i] += alpha * sum;
        }
    }
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
    const int M = *m_, N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    if (M <= 0 || N <= 0) return;

    /* Apply beta to C up front. Handles K==0 / alpha==0 paths in one
     * place — after this, the body computes only the matmul update. */
    int i, j;
    for (j = 0; j < N; ++j) {
        T *cj = &c[(size_t)j * ldc];
        if (beta == 0.0L)      for (i = 0; i < M; ++i) cj[i]  = 0.0L;
        else if (beta != 1.0L) for (i = 0; i < M; ++i) cj[i] *= beta;
    }
    if (alpha == 0.0L || K == 0) return;

    init_blocks();
    const int MC = g_mc, KC = g_kc, NC = g_nc;

    /*
     * Parallelize the NC loop: each thread takes a disjoint column
     * band of C. Per-thread packing scratch (Ap, Bp) — no sharing,
     * no atomics, no critical sections. Reduction order differs
     * across thread counts; accepted by the 10-ulp consistency
     * tolerance (see doc/parallel-blas-20260513.md §8).
     */
#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        T *Ap = aligned_alloc(64, (size_t)MC * KC * sizeof(T));
        T *Bp = aligned_alloc(64, (size_t)KC * NC * sizeof(T));
        if (Ap && Bp) {
            int jc, pc, ic;
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (jc = 0; jc < N; jc += NC) {
                const int jb = (N - jc < NC) ? (N - jc) : NC;
                for (pc = 0; pc < K; pc += KC) {
                    const int pb = (K - pc < KC) ? (K - pc) : KC;
                    pack_B(b, ldb, pc, jc, pb, jb, tb, Bp);
                    for (ic = 0; ic < M; ic += MC) {
                        const int ib = (M - ic < MC) ? (M - ic) : MC;
                        pack_A(a, lda, ic, pc, ib, pb, ta, Ap);
                        inner_kernel(ib, jb, pb, alpha, Ap, Bp,
                                     &c[(size_t)jc * ldc + ic], ldc);
                    }
                }
            }
        }
        free(Ap);
        free(Bp);
    }
}
