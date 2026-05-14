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
    }
}
