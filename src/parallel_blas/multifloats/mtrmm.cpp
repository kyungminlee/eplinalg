/*
 * mtrmm — multifloats real (double-double) triangular multiply.
 *
 * B := alpha · op(A) · B (SIDE='L') or B := alpha · B · op(A) (SIDE='R'),
 * where op(A) ∈ {A, Aᵀ}. Scalar-core diagonal blocks + mgemm trailing
 * update with coarse-N parallelism. No SIMD diagonal kernel (in
 * mtrsm that path was 5–13% of total; the SIMD gemm trailing
 * dominates and we get that for free via the mgemm_ call).
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

#define MTRMM_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}

int g_nb_trmm = 0;
int trmm_nb(void) {
    if (g_nb_trmm == 0) g_nb_trmm = env_int("MTRMM_NB", 64);
    return g_nb_trmm;
}

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};

inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (T x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }

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

/* ── SIDE = 'L' column-range cores ──────────────────────────────── */

inline void mtrmm_lln_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int k = M - 1; k >= 0; --k) {
            if (!dd_iszero(B_(k, j))) {
                T temp = alpha * B_(k, j);
                for (int i = M - 1; i > k; --i)
                    B_(i, j) = B_(i, j) + temp * A_(i, k);
                if (nounit) temp = temp * A_(k, k);
                B_(k, j) = temp;
            }
        }
    }
}

inline void mtrmm_lun_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int k = 0; k < M; ++k) {
            if (!dd_iszero(B_(k, j))) {
                T temp = alpha * B_(k, j);
                for (int i = 0; i < k; ++i)
                    B_(i, j) = B_(i, j) + temp * A_(i, k);
                if (nounit) temp = temp * A_(k, k);
                B_(k, j) = temp;
            }
        }
    }
}

inline void mtrmm_llt_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = B_(i, j);
            if (nounit) t = t * A_(i, i);
            for (int k = i + 1; k < M; ++k) t = t + A_(k, i) * B_(k, j);
            B_(i, j) = alpha * t;
        }
    }
}

inline void mtrmm_lut_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = B_(i, j);
            if (nounit) t = t * A_(i, i);
            for (int k = 0; k < i; ++k) t = t + A_(k, i) * B_(k, j);
            B_(i, j) = alpha * t;
        }
    }
}

/* ── SIDE = 'R' row-range cores ─────────────────────────────────── */

inline void mtrmm_rln_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = 0; j < N; ++j) {
        T t = alpha;
        if (nounit) t = t * A_(j, j);
        if (!dd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, j) = B_(i, j) * t;
        for (int k = j + 1; k < N; ++k) {
            if (!dd_iszero(A_(k, j))) {
                const T akj = alpha * A_(k, j);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = B_(i, j) + akj * B_(i, k);
            }
        }
    }
}

inline void mtrmm_run_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = N - 1; j >= 0; --j) {
        T t = alpha;
        if (nounit) t = t * A_(j, j);
        if (!dd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, j) = B_(i, j) * t;
        for (int k = 0; k < j; ++k) {
            if (!dd_iszero(A_(k, j))) {
                const T akj = alpha * A_(k, j);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = B_(i, j) + akj * B_(i, k);
            }
        }
    }
}

inline void mtrmm_rlt_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = N - 1; k >= 0; --k) {
        for (int j = k + 1; j < N; ++j) {
            if (!dd_iszero(A_(j, k))) {
                const T ajk = alpha * A_(j, k);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = B_(i, j) + ajk * B_(i, k);
            }
        }
        T t = alpha;
        if (nounit) t = t * A_(k, k);
        if (!dd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, k) = B_(i, k) * t;
    }
}

inline void mtrmm_rut_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < k; ++j) {
            if (!dd_iszero(A_(j, k))) {
                const T ajk = alpha * A_(j, k);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = B_(i, j) + ajk * B_(i, k);
            }
        }
        T t = alpha;
        if (nounit) t = t * A_(k, k);
        if (!dd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, k) = B_(i, k) * t;
    }
}

/* ── Standalone OMP wrappers (small-M fallback). */

#ifdef _OPENMP
#define MTRMM_OMP_WRAP_L(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (N >= MTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = static_cast<int>(static_cast<long long>(N) * tid / nt);            \
                int je  = static_cast<int>(static_cast<long long>(N) * (tid + 1) / nt);      \
                core(js, je, M, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit); }            \
    }
#define MTRMM_OMP_WRAP_R(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (M >= MTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int is  = static_cast<int>(static_cast<long long>(M) * tid / nt);            \
                int ie  = static_cast<int>(static_cast<long long>(M) * (tid + 1) / nt);      \
                core(is, ie, N, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else { core(0, M, N, alpha, a, lda, b, ldb, nounit); }            \
    }
#else
#define MTRMM_OMP_WRAP_L(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, N, M, alpha, a, lda, b, ldb, nounit);                       \
    }
#define MTRMM_OMP_WRAP_R(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, M, N, alpha, a, lda, b, ldb, nounit);                       \
    }
#endif

MTRMM_OMP_WRAP_L(mtrmm_lln, mtrmm_lln_core)
MTRMM_OMP_WRAP_L(mtrmm_lun, mtrmm_lun_core)
MTRMM_OMP_WRAP_L(mtrmm_llt, mtrmm_llt_core)
MTRMM_OMP_WRAP_L(mtrmm_lut, mtrmm_lut_core)
MTRMM_OMP_WRAP_R(mtrmm_rln, mtrmm_rln_core)
MTRMM_OMP_WRAP_R(mtrmm_run, mtrmm_run_core)
MTRMM_OMP_WRAP_R(mtrmm_rlt, mtrmm_rlt_core)
MTRMM_OMP_WRAP_R(mtrmm_rut, mtrmm_rut_core)

/* ── Blocked SIDE='L' ────────────────────────────────────────────── */

enum trmm_variant_L { LLN, LUN, LLT, LUT };

void blocked_chunk_L(trmm_variant_L V, int j_start, int j_end,
                     int M, int nb, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    const int my_N = j_end - j_start;
    if (my_N <= 0) return;

    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    T *B_chunk = &B_(0, j_start);

    if (V == LLN) {
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            mtrmm_lln_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            if (ic > 0) {
                mgemm_(NN, NN, &ib, &my_N, &ic, &alpha,
                       &A_(ic, 0), &lda,
                       B_chunk, &ldb, &one_dd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            ic -= nb;
        }
    } else if (V == LUN) {
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            mtrmm_lun_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int j0 = ic + ib;
                mgemm_(NN, NN, &ib, &my_N, &trailing, &alpha,
                       &A_(ic, j0), &lda,
                       &B_chunk[j0], &ldb, &one_dd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
        }
    } else if (V == LLT) {
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            mtrmm_llt_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int i0 = ic + ib;
                mgemm_(TN, NN, &ib, &my_N, &trailing, &alpha,
                       &A_(i0, ic), &lda,
                       &B_chunk[i0], &ldb, &one_dd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
        }
    } else { /* LUT */
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            mtrmm_lut_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            if (ic > 0) {
                mgemm_(TN, NN, &ib, &my_N, &ic, &alpha,
                       &A_(0, ic), &lda,
                       B_chunk, &ldb, &one_dd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            ic -= nb;
        }
    }
}

void blocked_dispatch_L(trmm_variant_L V, int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trmm_nb();
#ifdef _OPENMP
    if (N >= MTRMM_OMP_MIN && omp_get_max_threads() > 1) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nt  = omp_get_num_threads();
            int js  = static_cast<int>(static_cast<long long>(N) * tid / nt);
            int je  = static_cast<int>(static_cast<long long>(N) * (tid + 1) / nt);
            blocked_chunk_L(V, js, je, M, nb, alpha, a, lda, b, ldb, nounit);
        }
        return;
    }
#endif
    blocked_chunk_L(V, 0, N, M, nb, alpha, a, lda, b, ldb, nounit);
}

/* ── Blocked SIDE='R' ────────────────────────────────────────────── */

enum trmm_variant_R { RLN, RUN, RLT, RUT };

void blocked_chunk_R(trmm_variant_R V, int i_start, int i_end,
                     int N, int nb, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    const int my_M = i_end - i_start;
    if (my_M <= 0) return;
    T *B_chunk = &B_(i_start, 0);

    if (V == RLN) {
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            mtrmm_rln_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            const int trailing = N - (jc + jb);
            if (trailing > 0) {
                const int k0 = jc + jb;
                mgemm_(NN, NN, &my_M, &jb, &trailing, &alpha,
                       &B_chunk[static_cast<std::size_t>(k0) * ldb], &ldb,
                       &A_(k0, jc), &lda, &one_dd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
        }
    } else if (V == RUN) {
        int jc = ((N - 1) / nb) * nb;
        while (jc >= 0) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            mtrmm_run_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            if (jc > 0) {
                mgemm_(NN, NN, &my_M, &jb, &jc, &alpha,
                       B_chunk, &ldb,
                       &A_(0, jc), &lda, &one_dd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
            jc -= nb;
        }
    } else if (V == RLT) {
        int jc = ((N - 1) / nb) * nb;
        while (jc >= 0) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            mtrmm_rlt_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            if (jc > 0) {
                mgemm_(NN, TN, &my_M, &jb, &jc, &alpha,
                       B_chunk, &ldb,
                       &A_(jc, 0), &lda, &one_dd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
            jc -= nb;
        }
    } else { /* RUT */
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            mtrmm_rut_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            const int trailing = N - (jc + jb);
            if (trailing > 0) {
                const int k0 = jc + jb;
                mgemm_(NN, TN, &my_M, &jb, &trailing, &alpha,
                       &B_chunk[static_cast<std::size_t>(k0) * ldb], &ldb,
                       &A_(jc, k0), &lda, &one_dd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
        }
    }
}

void blocked_dispatch_R(trmm_variant_R V, int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trmm_nb();
#ifdef _OPENMP
    if (M >= MTRMM_OMP_MIN && omp_get_max_threads() > 1) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nt  = omp_get_num_threads();
            int is  = static_cast<int>(static_cast<long long>(M) * tid / nt);
            int ie  = static_cast<int>(static_cast<long long>(M) * (tid + 1) / nt);
            blocked_chunk_R(V, is, ie, N, nb, alpha, a, lda, b, ldb, nounit);
        }
        return;
    }
#endif
    blocked_chunk_R(V, 0, M, N, nb, alpha, a, lda, b, ldb, nounit);
}

} /* anonymous namespace */

extern "C" void mtrmm_(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    T *b, const int *ldb_,
    std::size_t side_len, std::size_t uplo_len, std::size_t transa_len, std::size_t diag_len)
{
    (void)side_len; (void)uplo_len; (void)transa_len; (void)diag_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, ldb = *ldb_;
    const T alpha = *alpha_;
    const char SIDE = up(side);
    const char UPLO = up(uplo);
    char TR         = up(transa);
    if (TR == 'C') TR = 'T';   /* real DD: conj-trans ≡ trans */
    const int nounit = (up(diag) != 'U');

    if (M == 0 || N == 0) return;

    if (dd_iszero(alpha)) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = zero_dd;
        return;
    }

    const int nb = trmm_nb();

    if (SIDE == 'L') {
        const int use_blocked = (M >= 2 * nb);
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_L(LLN, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_L(LUN, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_lun(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_L(LLT, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_L(LUT, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_lut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    } else {
        const int use_blocked = (N >= 2 * nb);
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_R(RLN, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_rln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_R(RUN, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_run(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_R(RLT, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_rlt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_R(RUT, M, N, alpha, a, lda, b, ldb, nounit);
                else             mtrmm_rut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    }
}

#undef A_
#undef B_
