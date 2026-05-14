/*
 * qgemm — kind16 (REAL(KIND=16) / __float128) GEMM overlay.
 *
 * Goto-style cache-blocked + packed, OpenMP-parallel across the
 * outer NC loop. Same structure as the kind10 egemm; only the
 * scalar type changes.
 *
 * __float128 has no hardware support on x86-64; every add/multiply
 * lowers to a libquadmath call (__addtf3, __multf3, ...). Each op
 * is hundreds of cycles. That makes the arithmetic strictly the
 * bottleneck — cache blocking and packing buy little, but OpenMP
 * across NC scales nearly linearly with cores.
 *
 * Fortran ABI:
 *   - subroutine name lowercased + trailing underscore: `qgemm_`
 *   - REAL(KIND=16) ↔ `__float128` (16-byte, 8-byte aligned in
 *     gfortran's COMPLEX/REAL slots; the C type is `__float128`
 *     from <quadmath.h>).
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef QBLAS_INLINE_SOFTFP
/* Inline __multf3/__addtf3 via vendored soft-fp template headers
 * instead of going through libgcc's external entry points. See
 * src/parallel_blas/kind16/qmath_inline.h. */
#include "qmath_inline.h"
#endif

typedef __float128 T;

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
    /* Same defaults as kind10 — autotune lands later. */
    g_mc = env_int("QBLAS_MC",  64);
    g_kc = env_int("QBLAS_KC", 128);
    g_nc = env_int("QBLAS_NC", 256);
}

static int trans_code(const char *p, size_t len) {
    (void)len;
    char c = (char)toupper((unsigned char)*p);
    return (c == 'C') ? 'T' : c;  /* real type: C == T */
}

/* ── Packers ──────────────────────────────────────────────────── */

/* Row-major Ap layout for inner-product (DDOT) kernel: see egemm.c. */
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
        for (p = 0; p < pb; ++p) {
            const T *src = &B[(size_t)(pc + p) * ldb + jc];
            for (j = 0; j < jb; ++j) Bp[(size_t)j * pb + p] = src[j];
        }
    }
}

/* ── Inner kernel ─────────────────────────────────────────────── */

/* Inner-product micro-kernel: Ap row-major, Bp col-major; one dot
 * product per (i,j) with accumulator in register.
 *
 * Under -DQBLAS_INLINE_SOFTFP the soft-float __multf3/__addtf3
 * bodies inline (see qmath_inline.h). Without the macro we go
 * through libgcc's external entry points via the `*`/`+` operators
 * — matches the historical behavior. */
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
            T sum = 0.0Q;
            for (p = 0; p < pb; ++p) {
#if defined(QBLAS_INLINE_SOFTFP)
                sum = qadd(sum, qmul(ai[p], bj[p]));
#elif defined(QBLAS_USE_FMAQ)
                /* Fused multiply-add: one libquadmath call (fmaq),
                 * one rounding. Replaces `sum += a*b` (two calls,
                 * two roundings via __multf3 + __addtf3). */
                sum = fmaq(ai[p], bj[p], sum);
#else
                sum += ai[p] * bj[p];
#endif
            }
#if defined(QBLAS_INLINE_SOFTFP)
            cj[i] = qadd(cj[i], qmul(alpha, sum));
#elif defined(QBLAS_USE_FMAQ)
            cj[i] = fmaq(alpha, sum, cj[i]);
#else
            cj[i] += alpha * sum;
#endif
        }
    }
}

/* ── Entry point ──────────────────────────────────────────────── */

void qgemm_(
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

    const T zero = 0.0Q, one = 1.0Q;
    int i, j;
    for (j = 0; j < N; ++j) {
        T *cj = &c[(size_t)j * ldc];
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
