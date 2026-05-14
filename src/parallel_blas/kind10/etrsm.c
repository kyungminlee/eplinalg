/*
 * etrsm — kind10 (REAL(KIND=10) / `long double`) triangular solve.
 *
 * Solves one of:
 *   op(A) · X = alpha · B          (SIDE='L')
 *   X · op(A) = alpha · B          (SIDE='R')
 *
 * where op(A) ∈ {A, Aᵀ, Aᴴ}; for real types Aᴴ ≡ Aᵀ. A is M×M (or
 * N×N) triangular (upper or lower; optionally unit-diagonal). B is
 * overwritten with the solution X.
 *
 * Stage 1 implementation: scalar unblocked, all 16 distinct algorithm
 * variants. Same rank-1-update loop structure as the upstream Netlib
 * reference DTRSM so the numerical behavior matches the migrated
 * archive to a tight tolerance.
 *
 * A blocked version (with `egemm` trailing updates) lands in a follow-
 * up commit. For unblocked the kernel is O(M²) per RHS column; the
 * blocked form gets the GEMM acceleration.
 *
 * Fortran ABI: same shape as the migrated entry point. Character args
 * have hidden trailing size_t length values appended; we ignore the
 * lengths and look at the first byte (upper-cased).
 */

#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Threshold below which OMP parallel-for on the column axis isn't
 * worth the parallel-region setup. */
#define ETRSM_OMP_N_MIN 32

typedef long double T;

/* Block size for the SIDE='L' blocked path. Env-tunable; default
 * picked empirically to match egemm's natural panel sizing. */
static int env_int(const char *name, int dflt) {
    const char *s = getenv(name);
    if (!s || !*s) return dflt;
    int v = atoi(s);
    return v > 0 ? v : dflt;
}
static int g_nb_trsm = 0;
static int trsm_nb(void) {
    if (g_nb_trsm == 0) g_nb_trsm = env_int("ETRSM_NB", 64);
    return g_nb_trsm;
}

/* Local egemm declaration (the overlay's own egemm_ — symbol is in
 * our static archive). We call it for trailing-matrix updates. */
extern void egemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const T *alpha,
    const T *a, const int *lda,
    const T *b, const int *ldb,
    const T *beta,
    T *c, const int *ldc,
    size_t transa_len, size_t transb_len);

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

/* Convenience accessors over column-major A and B. */
#define A_(i, j)  a[(size_t)(j) * lda + (i)]
#define B_(i, j)  b[(size_t)(j) * ldb + (i)]

/* ── SIDE = 'L': solve op(A) X = α B, A is M×M, B is M×N ───────── */

/* (L, L, N): forward substitution. */
static void trsm_lln(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
#ifdef _OPENMP
    #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < M; ++k) {
            if (B_(k, j) != 0.0L) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = k + 1; i < M; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

/* (L, U, N): back substitution. */
static void trsm_lun(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
#ifdef _OPENMP
    #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = M - 1; k >= 0; --k) {
            if (B_(k, j) != 0.0L) {
                if (nounit) B_(k, j) /= A_(k, k);
                const T bk = B_(k, j);
                for (int i = 0; i < k; ++i)
                    B_(i, j) -= bk * A_(i, k);
            }
        }
    }
}

/* (L, L, T): solve Aᵀ X = α B, A lower → Aᵀ upper, back sub on Aᵀ. */
static void trsm_llt(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    /* Reference: inner-product form, walks i from M-1 down to 0,
     * accumulating A[k,i]·B[k,j] for k > i (rows of A below diagonal). */
#ifdef _OPENMP
    #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        for (int i = M - 1; i >= 0; --i) {
            T t = alpha * B_(i, j);
            for (int k = i + 1; k < M; ++k) t -= A_(k, i) * B_(k, j);
            if (nounit) t /= A_(i, i);
            B_(i, j) = t;
        }
    }
}

/* (L, U, T): solve Aᵀ X = α B, A upper → Aᵀ lower, forward sub on Aᵀ. */
static void trsm_lut(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
#ifdef _OPENMP
    #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < M; ++i) {
            T t = alpha * B_(i, j);
            for (int k = 0; k < i; ++k) t -= A_(k, i) * B_(k, j);
            if (nounit) t /= A_(i, i);
            B_(i, j) = t;
        }
    }
}

/* ── Blocked SIDE='L' variants: alpha pre-scaling + egemm trailing
 * updates. The unblocked subkernels above are reused for each
 * `nb × nb` diagonal block. The flop savings come from doing the
 * trailing update via the optimized egemm rather than the
 * unblocked rank-1 update form. */

static void prescale_B(int M, int N, T alpha, T *b, int ldb)
{
    if (alpha == 1.0L) return;
    if (alpha == 0.0L) {
#ifdef _OPENMP
        #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = 0.0L;
        return;
    }
#ifdef _OPENMP
    #pragma omp parallel for if(N >= ETRSM_OMP_N_MIN) schedule(static)
#endif
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
}

/* (L, L, N) blocked: forward sub.
 *   For ic in 0, nb, 2nb, ...:
 *     B[ic:ic+ib, :] -= A[ic:ic+ib, 0:ic] · B[0:ic, :]   (egemm)
 *     solve A[ic:ic+ib, ic:ic+ib] · X = B[ic:ic+ib, :]   (unblocked)
 */
static void blocked_lln(int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
    prescale_B(M, N, alpha, b, ldb);
    const T m_one = -1.0L, one = 1.0L;
    const char NN[1] = {'N'};
    for (int ic = 0; ic < M; ic += nb) {
        const int ib = (M - ic < nb) ? (M - ic) : nb;
        if (ic > 0) {
            /* B[ic:ic+ib, :] -= A[ic:ic+ib, 0:ic] · B[0:ic, :] */
            egemm_(NN, NN, &ib, &N, &ic, &m_one,
                   &A_(ic, 0), &lda,
                   &B_(0, 0), &ldb, &one,
                   &B_(ic, 0), &ldb, 1, 1);
        }
        trsm_lln(ib, N, 1.0L, &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
    }
}

/* (L, U, N) blocked: back sub.
 *   For ic in M-nb..0 step -nb:
 *     B[ic:ic+ib, :] -= A[ic:ic+ib, ic+ib:M] · B[ic+ib:M, :]
 *     solve A[ic:ic+ib, ic:ic+ib] · X = B[ic:ic+ib, :]
 */
static void blocked_lun(int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
    prescale_B(M, N, alpha, b, ldb);
    const T m_one = -1.0L, one = 1.0L;
    const char NN[1] = {'N'};
    int ic = ((M - 1) / nb) * nb;     /* start of last block */
    while (ic >= 0) {
        const int ib = (M - ic < nb) ? (M - ic) : nb;
        const int trailing_rows = M - (ic + ib);
        if (trailing_rows > 0) {
            const int j0 = ic + ib;
            egemm_(NN, NN, &ib, &N, &trailing_rows, &m_one,
                   &A_(ic, j0), &lda,
                   &B_(j0, 0), &ldb, &one,
                   &B_(ic, 0), &ldb, 1, 1);
        }
        trsm_lun(ib, N, 1.0L, &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
        ic -= nb;
    }
}

/* (L, L, T) blocked: solve Aᵀ X = α B on lower-tri A.
 * Aᵀ is upper-tri; iterate top to bottom — back-sub-like but on Aᵀ.
 *   For ic in M-nb..0 step -nb:
 *     B[ic:ic+ib, :] -= Aᵀ[ic:ic+ib, ic+ib:M] · B[ic+ib:M, :]
 *                    = A[ic+ib:M, ic:ic+ib]ᵀ · B[ic+ib:M, :]
 *     solve Aᵀ[ic:ic+ib, ic:ic+ib] · X = B[ic:ic+ib, :]
 */
static void blocked_llt(int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
    prescale_B(M, N, alpha, b, ldb);
    const T m_one = -1.0L, one = 1.0L;
    const char TN[1] = {'T'};
    const char NN[1] = {'N'};
    int ic = ((M - 1) / nb) * nb;
    while (ic >= 0) {
        const int ib = (M - ic < nb) ? (M - ic) : nb;
        const int trailing_rows = M - (ic + ib);
        if (trailing_rows > 0) {
            const int i0 = ic + ib;
            /* GEMM with TRANSA='T': egemm('T','N', ib, N, trailing,
             *   -1, A[i0, ic] (k×m), lda, B[i0, 0], ldb, 1, B[ic,0], ldb)
             * Note the (k, m) shape because we're feeding A^T. */
            egemm_(TN, NN, &ib, &N, &trailing_rows, &m_one,
                   &A_(i0, ic), &lda,
                   &B_(i0, 0), &ldb, &one,
                   &B_(ic, 0), &ldb, 1, 1);
        }
        trsm_llt(ib, N, 1.0L, &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
        ic -= nb;
    }
}

/* (L, U, T) blocked: solve Aᵀ X = α B on upper-tri A.
 * Aᵀ is lower-tri; iterate top to bottom — forward-sub-like on Aᵀ.
 */
static void blocked_lut(int M, int N, T alpha,
                        const T *a, int lda, T *b, int ldb, int nounit)
{
    const int nb = trsm_nb();
    prescale_B(M, N, alpha, b, ldb);
    const T m_one = -1.0L, one = 1.0L;
    const char TN[1] = {'T'};
    const char NN[1] = {'N'};
    for (int ic = 0; ic < M; ic += nb) {
        const int ib = (M - ic < nb) ? (M - ic) : nb;
        if (ic > 0) {
            /* B[ic:ic+ib, :] -= Aᵀ[ic:ic+ib, 0:ic] · B[0:ic, :]
             *                 = A[0:ic, ic:ic+ib]ᵀ · B[0:ic, :] */
            egemm_(TN, NN, &ib, &N, &ic, &m_one,
                   &A_(0, ic), &lda,
                   &B_(0, 0), &ldb, &one,
                   &B_(ic, 0), &ldb, 1, 1);
        }
        trsm_lut(ib, N, 1.0L, &A_(ic, ic), lda, &B_(ic, 0), ldb, nounit);
    }
}

/* ── SIDE = 'R': solve X op(A) = α B, A is N×N, B is M×N ───────── */

/* (R, L, N): solve X · A = α B, A lower triangular.
 * Equivalent to back-substitution on columns of B from j = N-1 down. */
static void trsm_rln(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = N - 1; j >= 0; --j) {
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = j + 1; k < N; ++k) {
            if (A_(k, j) != 0.0L) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i)
                    B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = 1.0L / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) *= inv;
        }
    }
}

/* (R, U, N): solve X · A = α B, A upper triangular.
 * Forward-sub on columns of B from j = 0 up. */
static void trsm_run(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int j = 0; j < N; ++j) {
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, j) *= alpha;
        for (int k = 0; k < j; ++k) {
            if (A_(k, j) != 0.0L) {
                const T akj = A_(k, j);
                for (int i = 0; i < M; ++i)
                    B_(i, j) -= akj * B_(i, k);
            }
        }
        if (nounit) {
            const T inv = 1.0L / A_(j, j);
            for (int i = 0; i < M; ++i) B_(i, j) *= inv;
        }
    }
}

/* (R, L, T): solve X · Aᵀ = α B, A lower. */
static void trsm_rlt(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = 0; k < N; ++k) {
        if (nounit) {
            const T inv = 1.0L / A_(k, k);
            for (int i = 0; i < M; ++i) B_(i, k) *= inv;
        }
        for (int j = k + 1; j < N; ++j) {
            if (A_(j, k) != 0.0L) {
                const T ajk = A_(j, k);
                for (int i = 0; i < M; ++i)
                    B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, k) *= alpha;
    }
}

/* (R, U, T): solve X · Aᵀ = α B, A upper. */
static void trsm_rut(int M, int N, T alpha,
                     const T *a, int lda, T *b, int ldb, int nounit)
{
    for (int k = N - 1; k >= 0; --k) {
        if (nounit) {
            const T inv = 1.0L / A_(k, k);
            for (int i = 0; i < M; ++i) B_(i, k) *= inv;
        }
        for (int j = 0; j < k; ++j) {
            if (A_(j, k) != 0.0L) {
                const T ajk = A_(j, k);
                for (int i = 0; i < M; ++i)
                    B_(i, j) -= ajk * B_(i, k);
            }
        }
        if (alpha != 1.0L) for (int i = 0; i < M; ++i) B_(i, k) *= alpha;
    }
}

/* ── Entry point ──────────────────────────────────────────────── */

void etrsm_(
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

    /* alpha == 0 quick return: B becomes all zeros. */
    if (alpha == 0.0L) {
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < M; ++i) B_(i, j) = 0.0L;
        return;
    }

    if (SIDE == 'L') {
        /* Blocked path when M is large enough to amortize the egemm
         * call overhead. Threshold is twice the block size — below
         * that, blocking only adds noise. */
        const int use_blocked = (M >= 2 * trsm_nb());
        if (TR == 'N') {
            if (UPLO == 'L') {
                if (use_blocked) blocked_lln(M, N, alpha, a, lda, b, ldb, nounit);
                else             trsm_lln(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_lun(M, N, alpha, a, lda, b, ldb, nounit);
                else             trsm_lun(M, N, alpha, a, lda, b, ldb, nounit);
            }
        } else {
            if (UPLO == 'L') {
                if (use_blocked) blocked_llt(M, N, alpha, a, lda, b, ldb, nounit);
                else             trsm_llt(M, N, alpha, a, lda, b, ldb, nounit);
            } else {
                if (use_blocked) blocked_lut(M, N, alpha, a, lda, b, ldb, nounit);
                else             trsm_lut(M, N, alpha, a, lda, b, ldb, nounit);
            }
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
