/*
 * mgemm — multifloats real GEMM overlay (float64x2, double-double).
 *
 * C++ implementation using `multifloats::float64x2` with its overloaded
 * +, -, *, += operators. The kernel body is structurally identical to
 * the kind10/kind16 C paths; only the scalar type and arithmetic surface
 * change.
 *
 * Exported with extern "C" so the symbol is `mgemm_` (gfortran name
 * mangling) and the ABI is the POD `float64x2` (sizeof == 2*double),
 * matching gfortran's `type(real64x2)` sequence layout.
 *
 * Per-element cost: ~8 fp ops for an add, ~12 for a mul (Dekker /
 * Knuth EFTs). Arithmetic-bound; OMP scales near-linearly.
 */

#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef MBLAS_SIMD_DD
#include "mgemm_simd_kernel.h"
#endif

namespace mf = multifloats;
using T = mf::float64x2;

namespace {

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}

int g_mc = 0, g_kc = 0, g_nc = 0;
void init_blocks() {
    if (g_mc) return;
    g_mc = env_int("MBLAS_MC",  64);
    g_kc = env_int("MBLAS_KC", 128);
    g_nc = env_int("MBLAS_NC", 256);
}

int trans_code(const char *p, std::size_t /*len*/) {
    char c = static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
    return (c == 'C') ? 'T' : c;  /* real type: C == T */
}

void pack_A(const T * __restrict__ A, int lda,
            int ic, int pc, int ib, int pb,
            int ta, T * __restrict__ Ap)
{
    if (ta == 'N') {
        for (int p = 0; p < pb; ++p) {
            const T *src = &A[static_cast<std::size_t>(pc + p) * lda + ic];
            T *dst = &Ap[static_cast<std::size_t>(p) * ib];
            for (int i = 0; i < ib; ++i) dst[i] = src[i];
        }
    } else {
        for (int i = 0; i < ib; ++i) {
            const T *src = &A[static_cast<std::size_t>(ic + i) * lda + pc];
            for (int p = 0; p < pb; ++p)
                Ap[static_cast<std::size_t>(p) * ib + i] = src[p];
        }
    }
}

void pack_B(const T * __restrict__ B, int ldb,
            int pc, int jc, int pb, int jb,
            int tb, T * __restrict__ Bp)
{
    if (tb == 'N') {
        for (int j = 0; j < jb; ++j) {
            const T *src = &B[static_cast<std::size_t>(jc + j) * ldb + pc];
            T *dst = &Bp[static_cast<std::size_t>(j) * pb];
            for (int p = 0; p < pb; ++p) dst[p] = src[p];
        }
    } else {
        for (int p = 0; p < pb; ++p) {
            const T *src = &B[static_cast<std::size_t>(pc + p) * ldb + jc];
            for (int j = 0; j < jb; ++j)
                Bp[static_cast<std::size_t>(j) * pb + p] = src[j];
        }
    }
}

void inner_kernel(int ib, int jb, int pb, T alpha,
                  const T * __restrict__ Ap, const T * __restrict__ Bp,
                  T * __restrict__ C, int ldc)
{
    for (int j = 0; j < jb; ++j) {
        T *cj = &C[static_cast<std::size_t>(j) * ldc];
        const T *bj = &Bp[static_cast<std::size_t>(j) * pb];
        for (int p = 0; p < pb; ++p) {
            const T t = alpha * bj[p];
            const T *ap = &Ap[static_cast<std::size_t>(p) * ib];
            for (int i = 0; i < ib; ++i) cj[i] += t * ap[i];
        }
    }
}

#ifdef MBLAS_SIMD_DD

/*
 * SoA pack_B for SIMD path.
 * Layout: NR-column panels along the j-axis. Within each panel, the
 * (p, j_local) plane is stored as NR contiguous doubles per p in two
 * parallel arrays (hi, lo) — so loading 4 doubles into a ymm register
 * for one p iteration is a straight `vmovupd`.
 *
 *   Bp_hi[ panel*(pb*NR) + p*NR + j_local ] = op(B)[pc+p, jc+panel*NR+j_local].limbs[0]
 *   Bp_lo[ same offset ]                    = ... .limbs[1]
 *
 * Trailing panel (jb mod NR != 0) is zero-padded so the SIMD kernel
 * can always run the full NR-wide tile; the writeback masks padded
 * lanes.
 */
void pack_B_soa(const T * __restrict__ B, int ldb,
                int pc, int jc, int pb, int jb, int tb,
                double * __restrict__ Bp_hi, double * __restrict__ Bp_lo)
{
    constexpr int NR = simd_dd::NR;
    const int npanels = (jb + NR - 1) / NR;
    for (int panel = 0; panel < npanels; ++panel) {
        const int j0 = panel * NR;
        const int nr_eff = (jb - j0 < NR) ? (jb - j0) : NR;
        double *dst_hi = &Bp_hi[static_cast<std::size_t>(panel) * pb * NR];
        double *dst_lo = &Bp_lo[static_cast<std::size_t>(panel) * pb * NR];
        if (tb == 'N') {
            /* B[pc+p, jc+j0+c] = B[(jc+j0+c)*ldb + (pc+p)] */
            for (int c = 0; c < nr_eff; ++c) {
                const T *col = &B[static_cast<std::size_t>(jc + j0 + c) * ldb + pc];
                for (int p = 0; p < pb; ++p) {
                    dst_hi[p * NR + c] = col[p].limbs[0];
                    dst_lo[p * NR + c] = col[p].limbs[1];
                }
            }
            for (int c = nr_eff; c < NR; ++c)
                for (int p = 0; p < pb; ++p) {
                    dst_hi[p * NR + c] = 0.0;
                    dst_lo[p * NR + c] = 0.0;
                }
        } else {
            /* op(B)[p, j] = B[jc+j, pc+p] = B[(pc+p)*ldb + (jc+j)] */
            for (int p = 0; p < pb; ++p) {
                const T *row = &B[static_cast<std::size_t>(pc + p) * ldb + (jc + j0)];
                for (int c = 0; c < nr_eff; ++c) {
                    dst_hi[p * NR + c] = row[c].limbs[0];
                    dst_lo[p * NR + c] = row[c].limbs[1];
                }
                for (int c = nr_eff; c < NR; ++c) {
                    dst_hi[p * NR + c] = 0.0;
                    dst_lo[p * NR + c] = 0.0;
                }
            }
        }
    }
}

/*
 * SIMD writeback helper: alpha-scale a DD accumulator (acc_h, acc_l)
 * and merge into C[i, j0..j0+nr_eff-1].
 */
static inline __attribute__((always_inline)) void
simd_writeback(__m256d alpha_h, __m256d alpha_l,
               __m256d acc_h, __m256d acc_l,
               T *C_row_i, int ldc, int j0, int nr_eff)
{
    constexpr int NR = simd_dd::NR;
    __m256d ph, pl;
    simd_dd::dd_mul(alpha_h, alpha_l, acc_h, acc_l, ph, pl);
    alignas(32) double ph_a[NR], pl_a[NR];
    _mm256_store_pd(ph_a, ph);
    _mm256_store_pd(pl_a, pl);
    for (int j = 0; j < nr_eff; ++j) {
        T r;
        r.limbs[0] = ph_a[j];
        r.limbs[1] = pl_a[j];
        T &dst = C_row_i[static_cast<std::size_t>(j0 + j) * ldc];
        dst = dst + r;
    }
}

/*
 * SIMD inner kernel — MR=2 register tile.
 * 2 rows of A × NR=4 cols of B per micro-kernel call. The two rows'
 * dd_mul/dd_add chains run in parallel — same B load amortized across
 * both, gcc gets to interleave 8 FMAs per p iteration for ILP.
 *
 * Register budget: 4 acc (ach0/acl0/ach1/acl1) + 4 broadcasts
 * (ah0/al0/ah1/al1) + 2 B (bh/bl) + scratch ≈ 12 of 16 ymm regs.
 *
 * Trailing odd row (ib odd) handled by a separate MR=1 tail loop —
 * keeps the main loop branch-free.
 */
void inner_kernel_simd(int ib, int jb, int pb, T alpha,
                       const T * __restrict__ Ap,
                       const double * __restrict__ Bp_hi,
                       const double * __restrict__ Bp_lo,
                       T * __restrict__ C, int ldc)
{
    constexpr int NR = simd_dd::NR;
    const int j_panels = (jb + NR - 1) / NR;
    const __m256d alpha_h = _mm256_set1_pd(alpha.limbs[0]);
    const __m256d alpha_l = _mm256_set1_pd(alpha.limbs[1]);
    for (int jp = 0; jp < j_panels; ++jp) {
        const int j0 = jp * NR;
        const int nr_eff = (jb - j0 < NR) ? (jb - j0) : NR;
        const double *Bp_h_panel = &Bp_hi[static_cast<std::size_t>(jp) * pb * NR];
        const double *Bp_l_panel = &Bp_lo[static_cast<std::size_t>(jp) * pb * NR];
        int i = 0;
        /* MR=2 main loop */
        for (; i + 1 < ib; i += 2) {
            __m256d ach0 = _mm256_setzero_pd(), acl0 = _mm256_setzero_pd();
            __m256d ach1 = _mm256_setzero_pd(), acl1 = _mm256_setzero_pd();
            for (int p = 0; p < pb; ++p) {
                __m256d bh = _mm256_loadu_pd(&Bp_h_panel[p * NR]);
                __m256d bl = _mm256_loadu_pd(&Bp_l_panel[p * NR]);
                const T &a0 = Ap[static_cast<std::size_t>(p) * ib + i];
                const T &a1 = Ap[static_cast<std::size_t>(p) * ib + i + 1];
                __m256d ah0 = _mm256_set1_pd(a0.limbs[0]);
                __m256d al0 = _mm256_set1_pd(a0.limbs[1]);
                __m256d ah1 = _mm256_set1_pd(a1.limbs[0]);
                __m256d al1 = _mm256_set1_pd(a1.limbs[1]);
                __m256d rh0, rl0, rh1, rl1;
                simd_dd::dd_mul(ah0, al0, bh, bl, rh0, rl0);
                simd_dd::dd_mul(ah1, al1, bh, bl, rh1, rl1);
                simd_dd::dd_add(ach0, acl0, rh0, rl0, ach0, acl0);
                simd_dd::dd_add(ach1, acl1, rh1, rl1, ach1, acl1);
            }
            simd_writeback(alpha_h, alpha_l, ach0, acl0,
                           &C[i],     ldc, j0, nr_eff);
            simd_writeback(alpha_h, alpha_l, ach1, acl1,
                           &C[i + 1], ldc, j0, nr_eff);
        }
        /* MR=1 tail (only fires when ib is odd) */
        for (; i < ib; ++i) {
            __m256d acc_h = _mm256_setzero_pd();
            __m256d acc_l = _mm256_setzero_pd();
            for (int p = 0; p < pb; ++p) {
                const T &aval = Ap[static_cast<std::size_t>(p) * ib + i];
                __m256d ah = _mm256_set1_pd(aval.limbs[0]);
                __m256d al = _mm256_set1_pd(aval.limbs[1]);
                __m256d bh = _mm256_loadu_pd(&Bp_h_panel[p * NR]);
                __m256d bl = _mm256_loadu_pd(&Bp_l_panel[p * NR]);
                __m256d rh, rl;
                simd_dd::dd_mul(ah, al, bh, bl, rh, rl);
                simd_dd::dd_add(acc_h, acc_l, rh, rl, acc_h, acc_l);
            }
            simd_writeback(alpha_h, alpha_l, acc_h, acc_l,
                           &C[i], ldc, j0, nr_eff);
        }
    }
}

#endif /* MBLAS_SIMD_DD */

const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};

}  // namespace

extern "C" void mgemm_(
    const char *transa, const char *transb,
    const int *m_, const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    std::size_t transa_len, std::size_t transb_len)
{
    const int M = *m_, N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    if (M <= 0 || N <= 0) return;

    for (int j = 0; j < N; ++j) {
        T *cj = &c[static_cast<std::size_t>(j) * ldc];
        if (beta == zero_dd) {
            for (int i = 0; i < M; ++i) cj[i] = zero_dd;
        } else if (beta != one_dd) {
            for (int i = 0; i < M; ++i) cj[i] = cj[i] * beta;
        }
    }
    if (alpha == zero_dd || K == 0) return;

    init_blocks();
    const int MC = g_mc, KC = g_kc, NC = g_nc;

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        T *Ap = static_cast<T *>(std::aligned_alloc(
            64, static_cast<std::size_t>(MC) * KC * sizeof(T)));
#ifdef MBLAS_SIMD_DD
        /* SoA Bp: round NC up to NR boundary for the trailing panel. */
        constexpr int NR_simd = simd_dd::NR;
        const int NC_pad = ((NC + NR_simd - 1) / NR_simd) * NR_simd;
        double *Bp_hi = static_cast<double *>(std::aligned_alloc(
            64, static_cast<std::size_t>(KC) * NC_pad * sizeof(double)));
        double *Bp_lo = static_cast<double *>(std::aligned_alloc(
            64, static_cast<std::size_t>(KC) * NC_pad * sizeof(double)));
        if (Ap && Bp_hi && Bp_lo) {
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (int jc = 0; jc < N; jc += NC) {
                const int jb = (N - jc < NC) ? (N - jc) : NC;
                for (int pc = 0; pc < K; pc += KC) {
                    const int pb = (K - pc < KC) ? (K - pc) : KC;
                    pack_B_soa(b, ldb, pc, jc, pb, jb, tb, Bp_hi, Bp_lo);
                    for (int ic = 0; ic < M; ic += MC) {
                        const int ib = (M - ic < MC) ? (M - ic) : MC;
                        pack_A(a, lda, ic, pc, ib, pb, ta, Ap);
                        inner_kernel_simd(ib, jb, pb, alpha, Ap, Bp_hi, Bp_lo,
                                          &c[static_cast<std::size_t>(jc) * ldc + ic],
                                          ldc);
                    }
                }
            }
        }
        std::free(Ap);
        std::free(Bp_hi);
        std::free(Bp_lo);
#else
        T *Bp = static_cast<T *>(std::aligned_alloc(
            64, static_cast<std::size_t>(KC) * NC * sizeof(T)));
        if (Ap && Bp) {
#ifdef _OPENMP
            #pragma omp for schedule(static)
#endif
            for (int jc = 0; jc < N; jc += NC) {
                const int jb = (N - jc < NC) ? (N - jc) : NC;
                for (int pc = 0; pc < K; pc += KC) {
                    const int pb = (K - pc < KC) ? (K - pc) : KC;
                    pack_B(b, ldb, pc, jc, pb, jb, tb, Bp);
                    for (int ic = 0; ic < M; ic += MC) {
                        const int ib = (M - ic < MC) ? (M - ic) : MC;
                        pack_A(a, lda, ic, pc, ib, pb, ta, Ap);
                        inner_kernel(ib, jb, pb, alpha, Ap, Bp,
                                     &c[static_cast<std::size_t>(jc) * ldc + ic],
                                     ldc);
                    }
                }
            }
        }
        std::free(Ap);
        std::free(Bp);
#endif /* MBLAS_SIMD_DD */
    }
}
