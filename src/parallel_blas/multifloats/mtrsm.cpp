/*
 * mtrsm — multifloats real (double-double) triangular solve.
 *
 * Solves one of:
 *   op(A) · X = α · B          (SIDE='L')
 *   X · op(A) = α · B          (SIDE='R')
 *
 * where op(A) ∈ {A, Aᵀ, Aᴴ}; for real DD types Aᴴ ≡ Aᵀ. A is M×M
 * (or N×N) triangular (upper or lower; optionally unit-diagonal).
 * B is overwritten with the solution X.
 *
 * Implementation stages (this file lands stages 1+2 together, with
 * stage 3 — SIMD diagonal kernel — in a follow-up):
 *
 *   1. Scalar unblocked, all 16 distinct variants. C++ scalar code
 *      with multifloats overloaded operators inlined into the hot
 *      loop — already beats migrated heavily because the migrated
 *      mtrsm goes through gfortran elementals (one call per DD op).
 *
 *   2. Blocked SIDE='L' (4 variants) with mgemm trailing update +
 *      coarse-N parallelism (one outer omp parallel, threads
 *      partition columns of B). mgemm inside each thread runs with
 *      a 1-thread inner team due to OMP_NESTED=false default — so
 *      we get SIMD GEMM trailing-update for free.
 *
 *   3. SIMD 4-wide AVX2 diagonal kernel — separate follow-up.
 *
 * Fortran ABI: extern "C" symbol `mtrsm_`. Character args have
 * hidden trailing size_t lengths.
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

#define MTRSM_OMP_N_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}

int g_nb_trsm = 0;
int trsm_nb(void) {
    if (g_nb_trsm == 0) g_nb_trsm = env_int("MTRSM_NB", 64);
    return g_nb_trsm;
}

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};

inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (T x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }

/* mgemm extern — we call it for trailing-matrix updates. */
extern "C" void mgemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const T *alpha,
    const T *a, const int *lda,
    const T *b, const int *ldb,
    const T *beta,
    T *c, const int *ldc,
    std::size_t transa_len, std::size_t transb_len);

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]
#define B_(i, j)  b[static_cast<std::size_t>(j) * ldb + (i)]

/* ── Column-range "core" kernels: serial work over j ∈ [j_start, j_end). */

inline void mtrsm_lln_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * alpha;
        for (int k = 0; k < M; ++k) {
            if (!dd_iszero(B_(k, j))) {
                if (nounit) B_(k, j) = B_(k, j) / A_(k, k);
                const T bk = B_(k, j);
                for (int i = k + 1; i < M; ++i)
                    B_(i, j) = B_(i, j) - bk * A_(i, k);
            }
        }
    }
}

inline void mtrsm_lun_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * alpha;
        for (int k = M - 1; k >= 0; --k) {
            if (!dd_iszero(B_(k, j))) {
                if (nounit) B_(k, j) = B_(k, j) / A_(k, k);
                const T bk = B_(k, j);
                for (int i = 0; i < k; ++i)
                    B_(i, j) = B_(i, j) - bk * A_(i, k);
            }
        }
    }
}

inline void mtrsm_llt_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = alpha * B_(i, j);
            for (int k = i + 1; k < M; ++k) t = t - A_(k, i) * B_(k, j);
            if (nounit) t = t / A_(i, i);
            B_(i, j) = t;
        }
    }
}

inline void mtrsm_lut_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = alpha * B_(i, j);
            for (int k = 0; k < i; ++k) t = t - A_(k, i) * B_(k, j);
            if (nounit) t = t / A_(i, i);
            B_(i, j) = t;
        }
    }
}

/* ── SIDE = 'R': solve X op(A) = α B. Always scalar unblocked.
 * Same code paths as etrsm; SIDE='R' is much less common in LAPACK
 * so blocking it is deferred. */

inline void mtrsm_rln_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    /* For SIDE='R', the algorithm reorders columns of B (j loop is
     * outermost over A's columns), so j_start/j_end indicate which
     * B-rows this thread owns. (Actually for SIDE='R' the natural
     * partition is rows of B, not columns.) Handled by callers. */
    (void)j_start; (void)j_end;
    /* Serial implementation matches reference DTRSM 'R','L','N': */
    const int N = j_end;  /* placeholder — see caller */
    for (int j = N - 1; j >= 0; --j) {
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * alpha;
        for (int k = j + 1; k < N; ++k) {
            if (!dd_iszero(A_(k, j))) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i)
                    B_(i, j) = B_(i, j) - akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = one_dd / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * inv;
        }
    }
}

inline void mtrsm_run_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    (void)j_start;
    const int N = j_end;
    for (int j = 0; j < N; ++j) {
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * alpha;
        for (int k = 0; k < j; ++k) {
            if (!dd_iszero(A_(k, j))) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i)
                    B_(i, j) = B_(i, j) - akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = one_dd / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * inv;
        }
    }
}

inline void mtrsm_rlt_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    (void)j_start;
    const int N = j_end;
    for (int k = 0; k < N; ++k) {
        if (nounit) {
            const T inv = one_dd / A_(k, k);
            for (int i = 0; i < M; ++i) B_(i, k) = B_(i, k) * inv;
        }
        for (int j = k + 1; j < N; ++j) {
            if (!dd_iszero(A_(j, k))) {
                const T ajk = A_(j, k);
                for (int i = 0; i < M; ++i)
                    B_(i, j) = B_(i, j) - ajk * B_(i, k);
            }
        }
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, k) = B_(i, k) * alpha;
    }
}

inline void mtrsm_rut_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    (void)j_start;
    const int N = j_end;
    for (int k = N - 1; k >= 0; --k) {
        if (nounit) {
            const T inv = one_dd / A_(k, k);
            for (int i = 0; i < M; ++i) B_(i, k) = B_(i, k) * inv;
        }
        for (int j = 0; j < k; ++j) {
            if (!dd_iszero(A_(j, k))) {
                const T ajk = A_(j, k);
                for (int i = 0; i < M; ++i)
                    B_(i, j) = B_(i, j) - ajk * B_(i, k);
            }
        }
        if (!dd_isone(alpha)) for (int i = 0; i < M; ++i) B_(i, k) = B_(i, k) * alpha;
    }
}

/* Standalone unblocked SIDE='L' entries: wrap _core in own parallel region. */
#ifdef _OPENMP
#define MTRSM_OMP_WRAP(name, core)                                          \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit)               \
    {                                                                       \
        if (N >= MTRSM_OMP_N_MIN && omp_get_max_threads() > 1) {            \
            _Pragma("omp parallel")                                         \
            {                                                               \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = static_cast<int>((long long)N * tid / nt);        \
                int je  = static_cast<int>((long long)N * (tid + 1) / nt);  \
                core(js, je, M, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else {                                                            \
            core(0, N, M, alpha, a, lda, b, ldb, nounit);                   \
        }                                                                   \
    }
#else
#define MTRSM_OMP_WRAP(name, core)                                          \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit)               \
    {                                                                       \
        core(0, N, M, alpha, a, lda, b, ldb, nounit);                       \
    }
#endif

MTRSM_OMP_WRAP(mtrsm_lln, mtrsm_lln_core)
MTRSM_OMP_WRAP(mtrsm_lun, mtrsm_lun_core)
MTRSM_OMP_WRAP(mtrsm_llt, mtrsm_llt_core)
MTRSM_OMP_WRAP(mtrsm_lut, mtrsm_lut_core)

/* SIDE='R' kernels just take the full N — no column partitioning. */
void mtrsm_rln(int M, int N, T alpha,
               const T *a, int lda, T *b, int ldb, int nounit) {
    mtrsm_rln_core(0, N, M, alpha, a, lda, b, ldb, nounit);
}
void mtrsm_run(int M, int N, T alpha,
               const T *a, int lda, T *b, int ldb, int nounit) {
    mtrsm_run_core(0, N, M, alpha, a, lda, b, ldb, nounit);
}
void mtrsm_rlt(int M, int N, T alpha,
               const T *a, int lda, T *b, int ldb, int nounit) {
    mtrsm_rlt_core(0, N, M, alpha, a, lda, b, ldb, nounit);
}
void mtrsm_rut(int M, int N, T alpha,
               const T *a, int lda, T *b, int ldb, int nounit) {
    mtrsm_rut_core(0, N, M, alpha, a, lda, b, ldb, nounit);
}

/* ── Blocked SIDE='L' variants: coarse-grain parallelism across N.
 *
 * One outer omp parallel partitions columns of B across threads.
 * Each thread runs serial blocked-TRSM on its own column chunk:
 *   - mgemm trailing update (auto runs single-threaded internally due
 *     to OMP_NESTED=false → no inner team. The SIMD micro-kernel
 *     still runs at full SIMD width per call.)
 *   - scalar core diagonal solve on the thread's column range.
 */

inline void prescale_chunk(int j_start, int j_end, int M, T alpha,
                           T *b, int ldb)
{
    if (dd_isone(alpha)) return;
    if (dd_iszero(alpha)) {
        for (int j = j_start; j < j_end; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = zero_dd;
        return;
    }
    for (int j = j_start; j < j_end; ++j)
        for (int i = 0; i < M; ++i) B_(i, j) = B_(i, j) * alpha;
}

enum trsm_variant { LLN, LUN, LLT, LUT };

void blocked_chunk(trsm_variant V, int j_start, int j_end,
                   int M, int nb, T alpha,
                   const T *a, int lda, T *b, int ldb, int nounit)
{
    const int my_N = j_end - j_start;
    if (my_N <= 0) return;
    prescale_chunk(j_start, j_end, M, alpha, b, ldb);

    const T m_one = T{-1.0, 0.0};
    const T one   = T{ 1.0, 0.0};
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    T *B_chunk = &B_(0, j_start);

    if (V == LLN) {
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            if (ic > 0) {
                mgemm_(NN, NN, &ib, &my_N, &ic, &m_one,
                       &A_(ic, 0), &lda,
                       B_chunk, &ldb, &one,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            mtrsm_lln_core(j_start, j_end, ib, one_dd,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
        }
    } else if (V == LUN) {
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int j0 = ic + ib;
                mgemm_(NN, NN, &ib, &my_N, &trailing, &m_one,
                       &A_(ic, j0), &lda,
                       &B_chunk[j0], &ldb, &one,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            mtrsm_lun_core(j_start, j_end, ib, one_dd,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            ic -= nb;
        }
    } else if (V == LLT) {
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int i0 = ic + ib;
                mgemm_(TN, NN, &ib, &my_N, &trailing, &m_one,
                       &A_(i0, ic), &lda,
                       &B_chunk[i0], &ldb, &one,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            mtrsm_llt_core(j_start, j_end, ib, one_dd,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            ic -= nb;
        }
    } else { /* LUT */
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            if (ic > 0) {
                mgemm_(TN, NN, &ib, &my_N, &ic, &m_one,
                       &A_(0, ic), &lda,
                       B_chunk, &ldb, &one,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            mtrsm_lut_core(j_start, j_end, ib, one_dd,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
        }
    }
}

void blocked_dispatch(trsm_variant V, int M, int N, T alpha,
                      const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
#ifdef _OPENMP
    if (N >= MTRSM_OMP_N_MIN && omp_get_max_threads() > 1) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nt  = omp_get_num_threads();
            int js  = static_cast<int>((long long)N * tid / nt);
            int je  = static_cast<int>((long long)N * (tid + 1) / nt);
            blocked_chunk(V, js, je, M, nb, alpha, a, lda, b, ldb, nounit);
        }
        return;
    }
#endif
    blocked_chunk(V, 0, N, M, nb, alpha, a, lda, b, ldb, nounit);
}

void blocked_lln(int M, int N, T alpha, const T *a, int lda, T *b, int ldb, int nounit) {
    blocked_dispatch(LLN, M, N, alpha, a, lda, b, ldb, nounit);
}
void blocked_lun(int M, int N, T alpha, const T *a, int lda, T *b, int ldb, int nounit) {
    blocked_dispatch(LUN, M, N, alpha, a, lda, b, ldb, nounit);
}
void blocked_llt(int M, int N, T alpha, const T *a, int lda, T *b, int ldb, int nounit) {
    blocked_dispatch(LLT, M, N, alpha, a, lda, b, ldb, nounit);
}
void blocked_lut(int M, int N, T alpha, const T *a, int lda, T *b, int ldb, int nounit) {
    blocked_dispatch(LUT, M, N, alpha, a, lda, b, ldb, nounit);
}

}  // namespace

extern "C" void mtrsm_(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    T *b, const int *ldb_,
    std::size_t side_len, std::size_t uplo_len,
    std::size_t transa_len, std::size_t diag_len)
{
    (void)side_len; (void)uplo_len; (void)transa_len; (void)diag_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, ldb = *ldb_;
    const T alpha = *alpha_;
    const char SIDE = up(side);
    const char UPLO = up(uplo);
    char TR = up(transa);
    if (TR == 'C') TR = 'T';
    const int nounit = (up(diag) != 'U');

    if (M == 0 || N == 0) return;

    if (dd_iszero(alpha)) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = zero_dd;
        return;
    }

    if (SIDE == 'L') {
        const int use_blocked = (M >= 2 * trsm_nb());
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_lln(M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrsm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_lun(M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrsm_lun(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else {
            if (UPLO == 'L') {
                if (use_blocked) blocked_llt(M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrsm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_lut(M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrsm_lut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    } else {
        if (TR == 'N') {
            if (UPLO == 'L') mtrsm_rln(M, N, alpha, a, lda, b, ldb, nounit);
            else             mtrsm_run(M, N, alpha, a, lda, b, ldb, nounit);
        } else {
            if (UPLO == 'L') mtrsm_rlt(M, N, alpha, a, lda, b, ldb, nounit);
            else             mtrsm_rut(M, N, alpha, a, lda, b, ldb, nounit);
        }
    }
}

#undef A_
#undef B_
