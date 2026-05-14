/*
 * ygemm — kind10 complex GEMM overlay
 *   (COMPLEX(KIND=10), x86-64 _Complex long double, 2 × 80-bit
 *    extended packed into a 32-byte struct matching gfortran's
 *    complex(10) ABI).
 *
 * Same Goto-style cache-blocked + packed + OpenMP scaffold as egemm.
 * The packers absorb both ordinary transpose (T) and conjugate
 * transpose (C); after packing the inner kernel sees uniform
 * column-major Ap/Bp regardless of source orientation.
 *
 * 'T' (transpose, no conj) and 'C' (conjugate transpose) diverge
 * only in the conj() applied during pack. For TRANSA='C' we conj
 * each A element while packing; same for B. The kernel itself is
 * the same complex FMA loop.
 *
 * Block sizes shared with the real path via the same env knobs
 * (EBLAS_MC / KC / NC) — complex elements are 2× the size of real,
 * so the effective cache footprint is 2× larger; that's accepted at
 * this stage and revisited at autotune.
 *
 * Fortran ABI:
 *   - subroutine name lowercased + trailing underscore: `ygemm_`
 *   - scalars by pointer; complex scalar passed as a pointer to a
 *     pair of long doubles (real then imag)
 *   - character args followed by hidden trailing `size_t` lengths
 *   - COMPLEX(KIND=10) ↔ `_Complex long double` (32 bytes on x86-64)
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <complex.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef _Complex long double C10;

/* ── Block sizes (shared knobs with egemm) ────────────────────── */

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
    return (char)toupper((unsigned char)*p);
}

/* ── Packers ──────────────────────────────────────────────────── */

/*
 * Pack op(A)(ic:ic+ib, pc:pc+pb) into Ap, column-major (ib, pb).
 * After packing: Ap[p*ib + i] = op(A)[i, p].
 *
 *   ta == 'N':  op(A)[i,p] =        A[ic+i, pc+p]
 *   ta == 'T':  op(A)[i,p] =        A[pc+p, ic+i]
 *   ta == 'C':  op(A)[i,p] = conj(  A[pc+p, ic+i] )
 *
 * Conjugation collapses into pack — the kernel sees plain
 * multiplications afterwards.
 */
static void pack_A(const C10 *restrict A, int lda,
                   int ic, int pc, int ib, int pb,
                   int ta, C10 *restrict Ap)
{
    int i, p;
    if (ta == 'N') {
        for (p = 0; p < pb; ++p) {
            const C10 *src = &A[(size_t)(pc + p) * lda + ic];
            C10 *dst = &Ap[(size_t)p * ib];
            for (i = 0; i < ib; ++i) dst[i] = src[i];
        }
    } else if (ta == 'T') {
        for (i = 0; i < ib; ++i) {
            const C10 *src = &A[(size_t)(ic + i) * lda + pc];
            for (p = 0; p < pb; ++p) Ap[(size_t)p * ib + i] = src[p];
        }
    } else {  /* 'C' */
        for (i = 0; i < ib; ++i) {
            const C10 *src = &A[(size_t)(ic + i) * lda + pc];
            for (p = 0; p < pb; ++p) Ap[(size_t)p * ib + i] = conjl(src[p]);
        }
    }
}

static void pack_B(const C10 *restrict B, int ldb,
                   int pc, int jc, int pb, int jb,
                   int tb, C10 *restrict Bp)
{
    int p, j;
    if (tb == 'N') {
        for (j = 0; j < jb; ++j) {
            const C10 *src = &B[(size_t)(jc + j) * ldb + pc];
            C10 *dst = &Bp[(size_t)j * pb];
            for (p = 0; p < pb; ++p) dst[p] = src[p];
        }
    } else if (tb == 'T') {
        for (p = 0; p < pb; ++p) {
            const C10 *src = &B[(size_t)(pc + p) * ldb + jc];
            for (j = 0; j < jb; ++j) Bp[(size_t)j * pb + p] = src[j];
        }
    } else {  /* 'C' */
        for (p = 0; p < pb; ++p) {
            const C10 *src = &B[(size_t)(pc + p) * ldb + jc];
            for (j = 0; j < jb; ++j) Bp[(size_t)j * pb + p] = conjl(src[j]);
        }
    }
}

/* ── Inner kernel ─────────────────────────────────────────────── */

/*
 * C(ib × jb) += alpha * Ap[ib, pb] * Bp[pb, jb]
 * All column-major, contiguous after packing. Same shape as the
 * real-precision inner kernel; gcc emits a sequence of complex FMAs
 * (each ~4 long-double mults + 2 adds — ~50-cycle compute per op).
 */
static void inner_kernel(int ib, int jb, int pb, C10 alpha,
                         const C10 *restrict Ap, const C10 *restrict Bp,
                         C10 *restrict C, int ldc)
{
    int i, j, p;
    for (j = 0; j < jb; ++j) {
        C10 *cj = &C[(size_t)j * ldc];
        const C10 *bj = &Bp[(size_t)j * pb];
        for (p = 0; p < pb; ++p) {
            const C10 t = alpha * bj[p];
            const C10 *ap = &Ap[(size_t)p * ib];
            for (i = 0; i < ib; ++i) cj[i] += t * ap[i];
        }
    }
}

/* ── Entry point ──────────────────────────────────────────────── */

void ygemm_(
    const char *transa, const char *transb,
    const int *m_, const int *n_, const int *k_,
    const C10 *alpha_,
    const C10 *a, const int *lda_,
    const C10 *b, const int *ldb_,
    const C10 *beta_,
    C10 *c, const int *ldc_,
    size_t transa_len, size_t transb_len)
{
    const int M = *m_, N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const C10 alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    if (M <= 0 || N <= 0) return;

    /* C <- beta * C up front. Handles K==0 / alpha==0 paths. */
    int i, j;
    const C10 zero = 0.0L + 0.0iL;
    const C10 one  = 1.0L + 0.0iL;
    for (j = 0; j < N; ++j) {
        C10 *cj = &c[(size_t)j * ldc];
        if (beta == zero)      for (i = 0; i < M; ++i) cj[i]  = zero;
        else if (beta != one)  for (i = 0; i < M; ++i) cj[i] *= beta;
    }
    if (alpha == zero || K == 0) return;

    init_blocks();
    const int MC = g_mc, KC = g_kc, NC = g_nc;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        C10 *Ap = aligned_alloc(64, (size_t)MC * KC * sizeof(C10));
        C10 *Bp = aligned_alloc(64, (size_t)KC * NC * sizeof(C10));
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
