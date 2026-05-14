/*
 * xtrsm — kind16 complex (COMPLEX(KIND=16) / `__complex128`)
 *         triangular solve.
 *
 * Same scaffolding as qtrsm and the multifloats wtrsm scalar paths,
 * with TRANSA='C' (conjugate transpose) handled as a distinct case
 * from 'T'. No SIMD path — every __float128 op goes through
 * libgcc / libquadmath (the kind16 perf wall analyzed in
 * doc/parallel-blas-20260513.md §10).
 *
 * The xgemm call covers the trailing-matrix update for SIDE='L'.
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define XTRSM_OMP_N_MIN 32

typedef __complex128 T;

static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}

static int g_nb_trsm = 0;
static int trsm_nb(void) {
    if (g_nb_trsm == 0) g_nb_trsm = env_int("XTRSM_NB", 64);
    return g_nb_trsm;
}

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

static const T ZERO = 0.0Q + 0.0Qi;
static const T ONE  = 1.0Q + 0.0Qi;
static const T NEG_ONE = -1.0Q + 0.0Qi;

static inline T cconj(T a) { return conjq(a); }

extern void xgemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const T *alpha,
    const T *a, const int *lda,
    const T *b, const int *ldb,
    const T *beta,
    T *c, const int *ldc,
    size_t transa_len, size_t transb_len);

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]

static inline T A_op(const T *a, int lda, int row, int col, int conj_flag) {
    return conj_flag ? cconj(A_(row, col)) : A_(row, col);
}

/* ── Scalar column-range cores ───────────────────────────────── */

static inline void xtrsm_lln_core(int j_start, int j_end, int M, T alpha,
                                  const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < M; ++k) {
            if (B_(k, j) != ZERO) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = k + 1; i < M; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

static inline void xtrsm_lun_core(int j_start, int j_end, int M, T alpha,
                                  const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = M - 1; k >= 0; --k) {
            if (B_(k, j) != ZERO) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = 0; i < k; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

/* (L, L, T or C): inner-product form on op(A)ᵀ where op may conj. */
static inline void xtrsm_llTC_core(int j_start, int j_end, int M, T alpha,
                                   const T *a, int lda, T *b, int ldb,
                                   int nounit, int conj_flag)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = alpha * B_(i, j);
            for (int k = i + 1; k < M; ++k) t -= A_op(a, lda, k, i, conj_flag) * B_(k, j);
            if (nounit) t /= A_op(a, lda, i, i, conj_flag);
            B_(i, j) = t;
        }
    }
}

/* (L, U, T or C). */
static inline void xtrsm_luTC_core(int j_start, int j_end, int M, T alpha,
                                   const T *a, int lda, T *b, int ldb,
                                   int nounit, int conj_flag)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = alpha * B_(i, j);
            for (int k = 0; k < i; ++k) t -= A_op(a, lda, k, i, conj_flag) * B_(k, j);
            if (nounit) t /= A_op(a, lda, i, i, conj_flag);
            B_(i, j) = t;
        }
    }
}

/* ── SIDE = 'R' cores (scalar, no coarse-N). */

static void xtrsm_rln_core(int M, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = N - 1; j >= 0; --j) {
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = j + 1; k < N; ++k) {
            if (A_(k, j) != ZERO) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i) B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = ONE / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) *= inv;
        }
    }
}

static void xtrsm_run_core(int M, int N, T alpha,
                           const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = 0; j < N; ++j) {
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < j; ++k) {
            if (A_(k, j) != ZERO) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i) B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = ONE / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) *= inv;
        }
    }
}

static void xtrsm_rlTC_core(int M, int N, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int k = 0; k < N; ++k) {
        if (nounit) {
            const T inv = ONE / A_op(a, lda, k, k, conj_flag);
            for (int i = 0; i < M; ++i) B_(i, k) *= inv;
        }
        for (int j = k + 1; j < N; ++j) {
            const T ajk = A_op(a, lda, j, k, conj_flag);
            if (ajk != ZERO) {
                for (int i = 0; i < M; ++i) B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, k) *= alpha;
    }
}

static void xtrsm_ruTC_core(int M, int N, T alpha,
                            const T *a, int lda, T *b, int ldb,
                            int nounit, int conj_flag)
{
    for (int k = N - 1; k >= 0; --k) {
        if (nounit) {
            const T inv = ONE / A_op(a, lda, k, k, conj_flag);
            for (int i = 0; i < M; ++i) B_(i, k) *= inv;
        }
        for (int j = 0; j < k; ++j) {
            const T ajk = A_op(a, lda, j, k, conj_flag);
            if (ajk != ZERO) {
                for (int i = 0; i < M; ++i) B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != ONE) for (int i = 0; i < M; ++i) B_(i, k) *= alpha;
    }
}

/* Standalone OMP wrappers for SIDE='L' core variants. */
#ifdef _OPENMP
#define XTRSM_OMP_WRAP_LLN_LUN(name, core)                                  \
    static void name(int M, int N, T alpha,                                 \
                     const T *a, int lda, T *b, int ldb, int nounit)        \
    {                                                                       \
        if (N >= XTRSM_OMP_N_MIN && omp_get_max_threads() > 1) {            \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = (int)((long long)N * tid / nt);                   \
                int je  = (int)((long long)N * (tid + 1) / nt);             \
                core(js, je, M, alpha, a, lda, b, ldb, nounit);             \
            }                                                               \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit); }            \
    }
#define XTRSM_OMP_WRAP_TC(name, core, cflag)                                \
    static void name(int M, int N, T alpha,                                 \
                     const T *a, int lda, T *b, int ldb, int nounit)        \
    {                                                                       \
        if (N >= XTRSM_OMP_N_MIN && omp_get_max_threads() > 1) {            \
            _Pragma("omp parallel") {                                       \
                int tid = omp_get_thread_num();                             \
                int nt  = omp_get_num_threads();                            \
                int js  = (int)((long long)N * tid / nt);                   \
                int je  = (int)((long long)N * (tid + 1) / nt);             \
                core(js, je, M, alpha, a, lda, b, ldb, nounit, cflag);      \
            }                                                               \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit, cflag); }     \
    }
#else
#define XTRSM_OMP_WRAP_LLN_LUN(name, core)                                  \
    static void name(int M, int N, T alpha,                                 \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        core(0, N, M, alpha, a, lda, b, ldb, nounit);                       \
    }
#define XTRSM_OMP_WRAP_TC(name, core, cflag)                                \
    static void name(int M, int N, T alpha,                                 \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        core(0, N, M, alpha, a, lda, b, ldb, nounit, cflag);                \
    }
#endif

XTRSM_OMP_WRAP_LLN_LUN(xtrsm_lln, xtrsm_lln_core)
XTRSM_OMP_WRAP_LLN_LUN(xtrsm_lun, xtrsm_lun_core)
XTRSM_OMP_WRAP_TC(xtrsm_llt, xtrsm_llTC_core, 0)
XTRSM_OMP_WRAP_TC(xtrsm_lut, xtrsm_luTC_core, 0)
XTRSM_OMP_WRAP_TC(xtrsm_llc, xtrsm_llTC_core, 1)
XTRSM_OMP_WRAP_TC(xtrsm_luc, xtrsm_luTC_core, 1)

/* SIDE='R' direct entries. */
static void xtrsm_rln(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_rln_core(M, N, a_, a, lda, b, ldb, n); }
static void xtrsm_run(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_run_core(M, N, a_, a, lda, b, ldb, n); }
static void xtrsm_rlt(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_rlTC_core(M, N, a_, a, lda, b, ldb, n, 0); }
static void xtrsm_rut(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_ruTC_core(M, N, a_, a, lda, b, ldb, n, 0); }
static void xtrsm_rlc(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_rlTC_core(M, N, a_, a, lda, b, ldb, n, 1); }
static void xtrsm_ruc(int M, int N, T a_, const T *a, int lda, T *b, int ldb, int n) { xtrsm_ruTC_core(M, N, a_, a, lda, b, ldb, n, 1); }

/* ── Blocked SIDE='L' variants. */

static void prescale_chunk(int j_start, int j_end, int M, T alpha, T *b, int ldb)
{
    if (alpha == ONE) return;
    if (alpha == ZERO) {
        for (int j = j_start; j < j_end; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = ZERO;
        return;
    }
    for (int j = j_start; j < j_end; ++j)
        for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
}

enum xtrsm_variant { XLLN, XLUN, XLLT, XLUT, XLLC, XLUC };

static void blocked_chunk(enum xtrsm_variant V, int j_start, int j_end,
                          int M, int nb, T alpha,
                          const T *a, int lda, T *b, int ldb, int nounit)
{
    const int my_N = j_end - j_start;
    if (my_N <= 0) return;
    prescale_chunk(j_start, j_end, M, alpha, b, ldb);

    const char NN[1] = {'N'};
    const char TN[1] = {'T'};
    const char CN[1] = {'C'};
    T *B_chunk = &B_(0, j_start);

    if (V == XLLN) {
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            if (ic > 0) {
                xgemm_(NN, NN, &ib, &my_N, &ic, &NEG_ONE,
                       &A_(ic, 0), &lda,
                       B_chunk, &ldb, &ONE,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            xtrsm_lln_core(j_start, j_end, ib, ONE,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
        }
    } else if (V == XLUN) {
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int j0 = ic + ib;
                xgemm_(NN, NN, &ib, &my_N, &trailing, &NEG_ONE,
                       &A_(ic, j0), &lda,
                       &B_chunk[j0], &ldb, &ONE,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            xtrsm_lun_core(j_start, j_end, ib, ONE,
                           &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
            ic -= nb;
        }
    } else if (V == XLLT || V == XLLC) {
        const int conj_flag = (V == XLLC) ? 1 : 0;
        const char *trans_gemm = conj_flag ? CN : TN;
        int ic = ((M - 1) / nb) * nb;
        while (ic >= 0) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            const int trailing = M - (ic + ib);
            if (trailing > 0) {
                const int i0 = ic + ib;
                xgemm_(trans_gemm, NN, &ib, &my_N, &trailing, &NEG_ONE,
                       &A_(i0, ic), &lda,
                       &B_chunk[i0], &ldb, &ONE,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            xtrsm_llTC_core(j_start, j_end, ib, ONE,
                            &A_(ic, ic), lda, &B_(ic, 0), ldb,
                            nounit, conj_flag);
            ic -= nb;
        }
    } else { /* XLUT or XLUC */
        const int conj_flag = (V == XLUC) ? 1 : 0;
        const char *trans_gemm = conj_flag ? CN : TN;
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            if (ic > 0) {
                xgemm_(trans_gemm, NN, &ib, &my_N, &ic, &NEG_ONE,
                       &A_(0, ic), &lda,
                       B_chunk, &ldb, &ONE,
                       &B_chunk[ic], &ldb, 1, 1);
            }
            xtrsm_luTC_core(j_start, j_end, ib, ONE,
                            &A_(ic, ic), lda, &B_(ic, 0), ldb,
                            nounit, conj_flag);
        }
    }
}

static void blocked_dispatch(enum xtrsm_variant V, int M, int N, T alpha,
                             const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
#ifdef _OPENMP
    if (N >= XTRSM_OMP_N_MIN && omp_get_max_threads() > 1) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nt  = omp_get_num_threads();
            int js  = (int)((long long)N * tid / nt);
            int je  = (int)((long long)N * (tid + 1) / nt);
            blocked_chunk(V, js, je, M, nb, alpha, a, lda, b, ldb, nounit);
        }
        return;
    }
#endif
    blocked_chunk(V, 0, N, M, nb, alpha, a, lda, b, ldb, nounit);
}

void xtrsm_(
    const char *side, const char *uplo, const char *transa, const char *diag,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    T *b, const int *ldb_,
    size_t side_len, size_t uplo_len, size_t transa_len, size_t diag_len)
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

    if (alpha == ZERO) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = ZERO;
        return;
    }

    if (SIDE == 'L') {
        const int use_blocked = (M >= 2 * trsm_nb());
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch(XLLN, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch(XLUN, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_lun(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else if (TR == 'T') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch(XLLT, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch(XLUT, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_lut(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else { /* 'C' */
            if (UPLO == 'L') {
                if (use_blocked) blocked_dispatch(XLLC, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_llc(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_dispatch(XLUC, M, N, alpha, a, lda, b, ldb, nounit);
                else             xtrsm_luc(M, N, alpha, a, lda, b, ldb, nounit);
            }
        }
    } else {
        if (TR == 'N') {
            if (UPLO == 'L') xtrsm_rln(M, N, alpha, a, lda, b, ldb, nounit);
            else             xtrsm_run(M, N, alpha, a, lda, b, ldb, nounit);
        } else if (TR == 'T') {
            if (UPLO == 'L') xtrsm_rlt(M, N, alpha, a, lda, b, ldb, nounit);
            else             xtrsm_rut(M, N, alpha, a, lda, b, ldb, nounit);
        } else {
            if (UPLO == 'L') xtrsm_rlc(M, N, alpha, a, lda, b, ldb, nounit);
            else             xtrsm_ruc(M, N, alpha, a, lda, b, ldb, nounit);
        }
    }
}

#undef A_
#undef B_
