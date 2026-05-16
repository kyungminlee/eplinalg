/*
 * qtrsm — kind16 (REAL(KIND=16) / `__float128`) triangular solve.
 *
 * Solves one of:
 *   op(A) · X = alpha · B          (SIDE='L')
 *   X · op(A) = alpha · B          (SIDE='R')
 *
 * where op(A) ∈ {A, Aᵀ, Aᴴ}; for real types Aᴴ ≡ Aᵀ. A is M×M (or
 * N×N) triangular (upper or lower; optionally unit-diagonal). B is
 * overwritten with the solution X.
 *
 * Implementation: unblocked Netlib reference algorithm with OpenMP
 * coarse-grain parallelism. SIDE='L' partitions columns of B across
 * threads; SIDE='R' partitions rows of B (in both cases, the
 * partition axis carries no cross-thread dependence).
 *
 * No blocking / no qgemm trailing update: at kind16, every op lowers
 * to a libquadmath call, so blocking adds dispatch overhead without
 * accelerating the arithmetic. The reference algorithm matches
 * migrated DTRSM single-thread, and the OMP wrappers add the
 * parallel-scaling layer on top.
 */

#include <stddef.h>
#include <ctype.h>
#include <quadmath.h>
#ifdef _OPENMP
#include <omp.h>
#include "../common/blas_omp.h"
#endif

/* Threshold below which the parallel-region setup isn't worth it. */
#define QTRSM_OMP_MIN 32

typedef __float128 T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]

/* ── SIDE = 'L' column-range cores ────────────────────────────────
 * Each thread owns a column slice [j_start, j_end) of B and runs the
 * full Netlib algorithm on its slice. A is read-only and shared. */

static inline void trsm_lln_core(int j_start, int j_end, int M, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (alpha != 1.0Q) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < M; ++k) {
            if (B_(k, j) != 0.0Q) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = k + 1; i < M; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

static inline void trsm_lun_core(int j_start, int j_end, int M, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        if (alpha != 1.0Q) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = M - 1; k >= 0; --k) {
            if (B_(k, j) != 0.0Q) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = 0; i < k; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

static inline void trsm_llt_core(int j_start, int j_end, int M, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = alpha * B_(i, j);
            for (int k = i + 1; k < M; ++k) t -= A_(k, i) * B_(k, j);
            if (nounit) t /= A_(i, i);
            B_(i, j) = t;
        }
    }
}

static inline void trsm_lut_core(int j_start, int j_end, int M, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = j_start; j < j_end; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = alpha * B_(i, j);
            for (int k = 0; k < i; ++k) t -= A_(k, i) * B_(k, j);
            if (nounit) t /= A_(i, i);
            B_(i, j) = t;
        }
    }
}

/* ── SIDE = 'R' row-range cores ────────────────────────────────────
 * Each thread owns a row slice [i_start, i_end) of B. The algorithm
 * runs unchanged on its row strip — every B access is via B_(i, ·)
 * with i in the strip; A is read-only. */

static inline void trsm_rln_core(int i_start, int i_end, int N, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = N - 1; j >= 0; --j) {
        if (alpha != 1.0Q) for (int i = i_start; i < i_end; ++i) B_(i, j) *= alpha;
        for (int k = j + 1; k < N; ++k) {
            if (A_(k, j) != 0.0Q) {
                const T akj = A_(k, j);
                for (int i = i_start; i < i_end; ++i) B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = 1.0Q / A_(j, j);
            for (int i = i_start; i < i_end; ++i) B_(i, j) *= inv;
        }
    }
}

static inline void trsm_run_core(int i_start, int i_end, int N, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = 0; j < N; ++j) {
        if (alpha != 1.0Q) for (int i = i_start; i < i_end; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < j; ++k) {
            if (A_(k, j) != 0.0Q) {
                const T akj = A_(k, j);
                for (int i = i_start; i < i_end; ++i) B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = 1.0Q / A_(j, j);
            for (int i = i_start; i < i_end; ++i) B_(i, j) *= inv;
        }
    }
}

static inline void trsm_rlt_core(int i_start, int i_end, int N, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = 0; k < N; ++k) {
        if (nounit) {
            const T inv = 1.0Q / A_(k, k);
            for (int i = i_start; i < i_end; ++i) B_(i, k) *= inv;
        }
        for (int j = k + 1; j < N; ++j) {
            if (A_(j, k) != 0.0Q) {
                const T ajk = A_(j, k);
                for (int i = i_start; i < i_end; ++i) B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != 1.0Q) for (int i = i_start; i < i_end; ++i) B_(i, k) *= alpha;
    }
}

static inline void trsm_rut_core(int i_start, int i_end, int N, T alpha,
                                 const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = N - 1; k >= 0; --k) {
        if (nounit) {
            const T inv = 1.0Q / A_(k, k);
            for (int i = i_start; i < i_end; ++i) B_(i, k) *= inv;
        }
        for (int j = 0; j < k; ++j) {
            if (A_(j, k) != 0.0Q) {
                const T ajk = A_(j, k);
                for (int i = i_start; i < i_end; ++i) B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != 1.0Q) for (int i = i_start; i < i_end; ++i) B_(i, k) *= alpha;
    }
}

/* ── OMP wrappers: one parallel region per call, manual partition. ── */

#ifdef _OPENMP
#define QTRSM_OMP_WRAP_L(name, core)                                       \
    static void name(int M, int N, T alpha,                                \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        if (N >= QTRSM_OMP_MIN && blas_omp_max_threads() > 1) {             \
            _Pragma("omp parallel") {                                      \
                int tid = omp_get_thread_num();                            \
                int nt  = omp_get_num_threads();                           \
                int js  = (int)((long long)N * tid / nt);                  \
                int je  = (int)((long long)N * (tid + 1) / nt);            \
                core(js, je, M, alpha, a, lda, b, ldb, nounit);            \
            }                                                              \
        } else { core(0, N, M, alpha, a, lda, b, ldb, nounit); }           \
    }
#define QTRSM_OMP_WRAP_R(name, core)                                       \
    static void name(int M, int N, T alpha,                                \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        if (M >= QTRSM_OMP_MIN && blas_omp_max_threads() > 1) {             \
            _Pragma("omp parallel") {                                      \
                int tid = omp_get_thread_num();                            \
                int nt  = omp_get_num_threads();                           \
                int is  = (int)((long long)M * tid / nt);                  \
                int ie  = (int)((long long)M * (tid + 1) / nt);            \
                core(is, ie, N, alpha, a, lda, b, ldb, nounit);            \
            }                                                              \
        } else { core(0, M, N, alpha, a, lda, b, ldb, nounit); }           \
    }
#else
#define QTRSM_OMP_WRAP_L(name, core)                                       \
    static void name(int M, int N, T alpha,                                \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        core(0, N, M, alpha, a, lda, b, ldb, nounit);                      \
    }
#define QTRSM_OMP_WRAP_R(name, core)                                       \
    static void name(int M, int N, T alpha,                                \
                     const T *a, int lda, T *b, int ldb, int nounit) {     \
        core(0, M, N, alpha, a, lda, b, ldb, nounit);                      \
    }
#endif

QTRSM_OMP_WRAP_L(trsm_lln, trsm_lln_core)
QTRSM_OMP_WRAP_L(trsm_lun, trsm_lun_core)
QTRSM_OMP_WRAP_L(trsm_llt, trsm_llt_core)
QTRSM_OMP_WRAP_L(trsm_lut, trsm_lut_core)
QTRSM_OMP_WRAP_R(trsm_rln, trsm_rln_core)
QTRSM_OMP_WRAP_R(trsm_run, trsm_run_core)
QTRSM_OMP_WRAP_R(trsm_rlt, trsm_rlt_core)
QTRSM_OMP_WRAP_R(trsm_rut, trsm_rut_core)

/* ── Entry point ──────────────────────────────────────────────── */

void qtrsm_(
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
    const char SIDE   = up(side);
    const char UPLO   = up(uplo);
    char TR           = up(transa);
    if (TR == 'C') TR = 'T';   /* real type: conj-trans ≡ trans */
    const int nounit = (up(diag) != 'U');

    if (M == 0 || N == 0) return;

    if (alpha == 0.0Q) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = 0.0Q;
        return;
    }

    if (SIDE == 'L') {
        if (TR == 'N') {
            if (UPLO == 'L') trsm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            else             trsm_lun(M, N, alpha, a, lda, b, ldb, nounit);
        } else {
            if (UPLO == 'L') trsm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            else             trsm_lut(M, N, alpha, a, lda, b, ldb, nounit);
        }
    } else {
        if (TR == 'N') {
            if (UPLO == 'L') trsm_rln(M, N, alpha, a, lda, b, ldb, nounit);
            else             trsm_run(M, N, alpha, a, lda, b, ldb, nounit);
        } else {
            if (UPLO == 'L') trsm_rlt(M, N, alpha, a, lda, b, ldb, nounit);
            else             trsm_rut(M, N, alpha, a, lda, b, ldb, nounit);
        }
    }
}

#undef A_
#undef B_
