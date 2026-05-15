/*
 * wtrmm — multifloats complex (complex64x2) triangular multiply.
 *
 * Direct port of mtrmm.cpp to complex DD. Scalar diagonal cores +
 * wgemm trailing update + coarse-grain OMP. TRANSA='C' is conjugate
 * transpose (distinct from 'T'). No SIMD diagonal kernel — wgemm
 * itself is SIMD so the diagonal stays scalar (matches wtrsm's
 * structure where SIMD diagonal was a small fraction of total time).
 */

#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {

#define WTRMM_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}

int g_nb_trmm = 0;
int trmm_nb(void) {
    if (g_nb_trmm == 0) g_nb_trmm = env_int("WTRMM_NB", 64);
    return g_nb_trmm;
}

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_cdd{ R{0.0, 0.0}, R{0.0, 0.0} };
const T one_cdd { R{1.0, 0.0}, R{0.0, 0.0} };

inline bool cdd_iszero(const T &x) {
    return x.re.limbs[0] == 0.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}
inline bool cdd_isone(const T &x) {
    return x.re.limbs[0] == 1.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}

inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
inline T cconj(T const &a) { return T{ a.re, -a.im }; }

extern "C" void wgemm_(
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

inline T A_op(const T *a, int lda, int row, int col, int conj_flag) {
    return conj_flag ? cconj(A_(row, col)) : A_(row, col);
}

/* ── SIDE = 'L' column-range cores ──────────────────────────────── */

inline void wtrmm_lln_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int k = M - 1; k >= 0; --k) {
            if (!cdd_iszero(B_(k, j))) {
                T temp = cmul(alpha, B_(k, j));
                for (int i = M - 1; i > k; --i)
                    B_(i, j) = cadd(B_(i, j), cmul(temp, A_(i, k)));
                if (nounit) temp = cmul(temp, A_(k, k));
                B_(k, j) = temp;
            }
        }
    }
}

inline void wtrmm_lun_core(int j_start, int j_end, int M, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int k = 0; k < M; ++k) {
            if (!cdd_iszero(B_(k, j))) {
                T temp = cmul(alpha, B_(k, j));
                for (int i = 0; i < k; ++i)
                    B_(i, j) = cadd(B_(i, j), cmul(temp, A_(i, k)));
                if (nounit) temp = cmul(temp, A_(k, k));
                B_(k, j) = temp;
            }
        }
    }
}

inline void wtrmm_llTC_core(int j_start, int j_end, int M, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = B_(i, j);
            if (nounit) t = cmul(t, A_op(a, lda, i, i, conj_flag));
            for (int k = i + 1; k < M; ++k)
                t = cadd(t, cmul(A_op(a, lda, k, i, conj_flag), B_(k, j)));
            B_(i, j) = cmul(alpha, t);
        }
    }
}

inline void wtrmm_luTC_core(int j_start, int j_end, int M, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = B_(i, j);
            if (nounit) t = cmul(t, A_op(a, lda, i, i, conj_flag));
            for (int k = 0; k < i; ++k)
                t = cadd(t, cmul(A_op(a, lda, k, i, conj_flag), B_(k, j)));
            B_(i, j) = cmul(alpha, t);
        }
    }
}

/* ── SIDE = 'R' row-range cores ─────────────────────────────────── */

inline void wtrmm_rln_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = 0; j < N; ++j) {
        T t = alpha;
        if (nounit) t = cmul(t, A_(j, j));
        if (!cdd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, j) = cmul(B_(i, j), t);
        for (int k = j + 1; k < N; ++k) {
            if (!cdd_iszero(A_(k, j))) {
                const T akj = cmul(alpha, A_(k, j));
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = cadd(B_(i, j), cmul(akj, B_(i, k)));
            }
        }
    }
}

inline void wtrmm_run_core(int i_start, int i_end, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = N - 1; j >= 0; --j) {
        T t = alpha;
        if (nounit) t = cmul(t, A_(j, j));
        if (!cdd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, j) = cmul(B_(i, j), t);
        for (int k = 0; k < j; ++k) {
            if (!cdd_iszero(A_(k, j))) {
                const T akj = cmul(alpha, A_(k, j));
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = cadd(B_(i, j), cmul(akj, B_(i, k)));
            }
        }
    }
}

inline void wtrmm_rlTC_core(int i_start, int i_end, int N, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int k = N - 1; k >= 0; --k) {
        for (int j = k + 1; j < N; ++j) {
            const T ajk = A_op(a, lda, j, k, conj_flag);
            if (!cdd_iszero(ajk)) {
                const T scaled = cmul(alpha, ajk);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = cadd(B_(i, j), cmul(scaled, B_(i, k)));
            }
        }
        T t = alpha;
        if (nounit) t = cmul(t, A_op(a, lda, k, k, conj_flag));
        if (!cdd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, k) = cmul(B_(i, k), t);
    }
}

inline void wtrmm_ruTC_core(int i_start, int i_end, int N, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int k = 0; k < N; ++k) {
        for (int j = 0; j < k; ++j) {
            const T ajk = A_op(a, lda, j, k, conj_flag);
            if (!cdd_iszero(ajk)) {
                const T scaled = cmul(alpha, ajk);
                for (int i = i_start; i < i_end; ++i)
                    B_(i, j) = cadd(B_(i, j), cmul(scaled, B_(i, k)));
            }
        }
        T t = alpha;
        if (nounit) t = cmul(t, A_op(a, lda, k, k, conj_flag));
        if (!cdd_isone(t))
            for (int i = i_start; i < i_end; ++i) B_(i, k) = cmul(B_(i, k), t);
    }
}

/* ── Standalone OMP wrappers (small-M/N fallback). */

#ifdef _OPENMP
#define WTRMM_OMP_WRAP_L(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (N >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = static_cast<int>(static_cast<long long>(N) * tid / nt);            \
                int je  = static_cast<int>(static_cast<long long>(N) * (tid + 1) / nt);      \
                core(js, je, M, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit); }            \
    }
#define WTRMM_OMP_WRAP_L_TC(name, core, cflag)                              \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (N >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = static_cast<int>(static_cast<long long>(N) * tid / nt);            \
                int je  = static_cast<int>(static_cast<long long>(N) * (tid + 1) / nt);      \
                core(js, je, M, alpha, a, lda, b, ldb, nounit, cflag);      \
            }                                                               \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit, cflag); }     \
    }
#define WTRMM_OMP_WRAP_R(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (M >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int is  = static_cast<int>(static_cast<long long>(M) * tid / nt);            \
                int ie  = static_cast<int>(static_cast<long long>(M) * (tid + 1) / nt);      \
                core(is, ie, N, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else { core(0, M, N, alpha, a, lda, b, ldb, nounit); }            \
    }
#define WTRMM_OMP_WRAP_R_TC(name, core, cflag)                              \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        if (M >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {              \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int is  = static_cast<int>(static_cast<long long>(M) * tid / nt);            \
                int ie  = static_cast<int>(static_cast<long long>(M) * (tid + 1) / nt);      \
                core(is, ie, N, alpha, a, lda, b, ldb, nounit, cflag);      \
            }                                                               \
        } else { core(0, M, N, alpha, a, lda, b, ldb, nounit, cflag); }     \
    }
#else
#define WTRMM_OMP_WRAP_L(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, N, M, alpha, a, lda, b, ldb, nounit);                       \
    }
#define WTRMM_OMP_WRAP_L_TC(name, core, cflag)                              \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, N, M, alpha, a, lda, b, ldb, nounit, cflag);                \
    }
#define WTRMM_OMP_WRAP_R(name, core)                                       \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, M, N, alpha, a, lda, b, ldb, nounit);                       \
    }
#define WTRMM_OMP_WRAP_R_TC(name, core, cflag)                              \
    void name(int M, int N, T alpha,                                        \
              const T *a, int lda, T *b, int ldb, int nounit) {             \
        core(0, M, N, alpha, a, lda, b, ldb, nounit, cflag);                \
    }
#endif

WTRMM_OMP_WRAP_L   (wtrmm_lln, wtrmm_lln_core)
WTRMM_OMP_WRAP_L   (wtrmm_lun, wtrmm_lun_core)
WTRMM_OMP_WRAP_L_TC(wtrmm_llt, wtrmm_llTC_core, 0)
WTRMM_OMP_WRAP_L_TC(wtrmm_lut, wtrmm_luTC_core, 0)
WTRMM_OMP_WRAP_L_TC(wtrmm_llc, wtrmm_llTC_core, 1)
WTRMM_OMP_WRAP_L_TC(wtrmm_luc, wtrmm_luTC_core, 1)
WTRMM_OMP_WRAP_R   (wtrmm_rln, wtrmm_rln_core)
WTRMM_OMP_WRAP_R   (wtrmm_run, wtrmm_run_core)
WTRMM_OMP_WRAP_R_TC(wtrmm_rlt, wtrmm_rlTC_core, 0)
WTRMM_OMP_WRAP_R_TC(wtrmm_rut, wtrmm_ruTC_core, 0)
WTRMM_OMP_WRAP_R_TC(wtrmm_rlc, wtrmm_rlTC_core, 1)
WTRMM_OMP_WRAP_R_TC(wtrmm_ruc, wtrmm_ruTC_core, 1)

/* ── Blocked SIDE='L' ────────────────────────────────────────────── */

enum trmm_variant_L { WLLN, WLUN, WLLT, WLUT, WLLC, WLUC };

void blocked_chunk_L(trmm_variant_L V, int j_start, int j_end,
                     int M, int nb, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    const int my_N = j_end - j_start;
    if (my_N <= 0) return;

    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    const char CN[1] = {'C'};
    T *B_chunk = &B_(0, j_start);

    if (V == WLLN) {
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            wtrmm_lln_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            if (ic > 0) {
                wgemm_(NN, NN, &ib, &my_N, &ic, &alpha,
                       &A_(ic, 0), &lda,
                       B_chunk, &ldb, &one_cdd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            ic -= nb;
        }
    } else if (V == WLUN) {
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            wtrmm_lun_core(j_start, j_end, ib, alpha,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int j0 = ic + ib;
                wgemm_(NN, NN, &ib, &my_N, &trailing, &alpha,
                       &A_(ic, j0), &lda,
                       &B_chunk[j0], &ldb, &one_cdd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
        }
    } else if (V == WLLT || V == WLLC) {
        const int conj_flag = (V == WLLC) ? 1 : 0;
        const char *gemm_trans = conj_flag ? CN : TN;
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            wtrmm_llTC_core(j_start, j_end, ib, alpha,
                            &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit, conj_flag);
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int i0 = ic + ib;
                wgemm_(gemm_trans, NN, &ib, &my_N, &trailing, &alpha,
                       &A_(i0, ic), &lda,
                       &B_chunk[i0], &ldb, &one_cdd,
                       &B_chunk[ic], &ldb, 1, 1);
            }
        }
    } else { /* WLUT or WLUC */
        const int conj_flag = (V == WLUC) ? 1 : 0;
        const char *gemm_trans = conj_flag ? CN : TN;
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            wtrmm_luTC_core(j_start, j_end, ib, alpha,
                            &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit, conj_flag);
            if (ic > 0) {
                wgemm_(gemm_trans, NN, &ib, &my_N, &ic, &alpha,
                       &A_(0, ic), &lda,
                       B_chunk, &ldb, &one_cdd,
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
    if (N >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {
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

enum trmm_variant_R { WRLN, WRUN, WRLT, WRUT, WRLC, WRUC };

void blocked_chunk_R(trmm_variant_R V, int i_start, int i_end,
                     int N, int nb, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    const char CN[1] = {'C'};
    const int my_M = i_end - i_start;
    if (my_M <= 0) return;
    T *B_chunk = &B_(i_start, 0);

    if (V == WRLN) {
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            wtrmm_rln_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            const int trailing = N - (jc + jb);
            if (trailing > 0) {
                const int k0 = jc + jb;
                wgemm_(NN, NN, &my_M, &jb, &trailing, &alpha,
                       &B_chunk[static_cast<std::size_t>(k0) * ldb], &ldb,
                       &A_(k0, jc), &lda, &one_cdd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
        }
    } else if (V == WRUN) {
        int jc = ((N - 1) / nb) * nb;
        while (jc >= 0) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            wtrmm_run_core(i_start, i_end, jb, alpha,
                           &A_(jc, jc), lda, &B_(0, jc), ldb, nounit);
            if (jc > 0) {
                wgemm_(NN, NN, &my_M, &jb, &jc, &alpha,
                       B_chunk, &ldb,
                       &A_(0, jc), &lda, &one_cdd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
            jc -= nb;
        }
    } else if (V == WRLT || V == WRLC) {
        const int conj_flag = (V == WRLC) ? 1 : 0;
        const char *gemm_trans = conj_flag ? CN : TN;
        int jc = ((N - 1) / nb) * nb;
        while (jc >= 0) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            wtrmm_rlTC_core(i_start, i_end, jb, alpha,
                            &A_(jc, jc), lda, &B_(0, jc), ldb, nounit, conj_flag);
            if (jc > 0) {
                wgemm_(NN, gemm_trans, &my_M, &jb, &jc, &alpha,
                       B_chunk, &ldb,
                       &A_(jc, 0), &lda, &one_cdd,
                       &B_chunk[static_cast<std::size_t>(jc) * ldb], &ldb, 1, 1);
            }
            jc -= nb;
        }
    } else { /* WRUT or WRUC */
        const int conj_flag = (V == WRUC) ? 1 : 0;
        const char *gemm_trans = conj_flag ? CN : TN;
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            wtrmm_ruTC_core(i_start, i_end, jb, alpha,
                            &A_(jc, jc), lda, &B_(0, jc), ldb, nounit, conj_flag);
            const int trailing = N - (jc + jb);
            if (trailing > 0) {
                const int k0 = jc + jb;
                wgemm_(NN, gemm_trans, &my_M, &jb, &trailing, &alpha,
                       &B_chunk[static_cast<std::size_t>(k0) * ldb], &ldb,
                       &A_(jc, k0), &lda, &one_cdd,
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
    if (M >= WTRMM_OMP_MIN && omp_get_max_threads() > 1) {
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

extern "C" void wtrmm_(
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
    const char TR = up(transa);
    const int nounit = (up(diag) != 'U');

    if (M == 0 || N == 0) return;

    if (cdd_iszero(alpha)) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = zero_cdd;
        return;
    }

    const int nb = trmm_nb();

    if (SIDE == 'L') {
        const int use_blocked = (M >= 2 * nb);
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_L(WLLN, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_L(WLUN, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_lun(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else if (TR == 'T') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_L(WLLT, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_L(WLUT, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_lut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else { /* 'C' */
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_L(WLLC, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_llc(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_L(WLUC, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_luc(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    } else {
        const int use_blocked = (N >= 2 * nb);
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_R(WRLN, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_rln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_R(WRUN, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_run(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else if (TR == 'T') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_R(WRLT, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_rlt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_R(WRUT, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_rut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch_R(WRLC, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_rlc(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch_R(WRUC, M, N, alpha, a, lda, b, ldb, nounit);
                else             wtrmm_ruc(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    }
}

#undef A_
#undef B_
