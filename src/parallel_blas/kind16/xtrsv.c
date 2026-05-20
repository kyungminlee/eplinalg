/*
 * xtrsv — kind16 complex triangular solve.
 */

#include <stddef.h>
#include <stdlib.h>
#include <ctype.h>
#include <quadmath.h>

typedef __complex128 T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

#define XTRSV_BLOCKED_NB_DEFAULT 64

static int xtrsv_blocked_nb(void);

/* Forward declaration — body below. The dispatch in xtrsv_ routes
 * large stride-1 calls into the blocked path; the blocked path
 * recursively calls xtrsv_ for each NB-sized diagonal block, but
 * those sub-calls have N < 2·NB so they fall through to the serial
 * body and there is no infinite recursion. */
void xtrsv_blocked_(
    const char *uplo, const char *trans, const char *diag,
    const int *n_,
    const T *restrict a, const int *lda_,
    T *restrict x, const int *incx_,
    size_t uplo_len, size_t trans_len, size_t diag_len);

void xtrsv_(
    const char *uplo, const char *trans, const char *diag,
    const int *n_,
    const T *restrict a, const int *lda_,
    T *restrict x, const int *incx_,
    size_t uplo_len, size_t trans_len, size_t diag_len)
{
    (void)uplo_len; (void)trans_len; (void)diag_len;
    const int N = *n_;
    const int lda = *lda_, incx = *incx_;
    const char UPLO = up(uplo);
    const char TR   = up(trans);
    const char DIAG = up(diag);
    const int nounit = (DIAG != 'U');

    if (N == 0) return;

    /* Block-parallel dispatch: stride-1 calls above the blocking
     * threshold route through the LAPACK-blocked path, which uses
     * parallel xgemv_ for trailing updates. At OMP=4 this delivers
     * ~3.3-3.7× of the unblocked algorithm at N=1024 across all
     * uplo/trans/diag combos (kind16 bench 2026-05-20). At OMP=1
     * overhead is within noise (±5% vs unblocked). */
    if (incx == 1 && N >= 2 * xtrsv_blocked_nb()) {
        xtrsv_blocked_(uplo, trans, diag, n_, a, lda_, x, incx_,
                       uplo_len, trans_len, diag_len);
        return;
    }

    const T zero = 0.0Q + 0.0Qi;

    if (incx == 1) {
        if (TR == 'N') {
            if (UPLO == 'L') {
                for (int i = 0; i < N; ++i) {
                    if (x[i] != zero) {
                        if (nounit) x[i] /= A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = i + 1; k < N; ++k) x[k] -= xi * ai[k];
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    if (x[i] != zero) {
                        if (nounit) x[i] /= A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = 0; k < i; ++k) x[k] -= xi * ai[k];
                    }
                }
            }
        } else {
            const int conj_a = (TR == 'C');
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    if (conj_a) {
                        for (int k = i + 1; k < N; ++k) t -= conjq(ai[k]) * x[k];
                        if (nounit) t /= conjq(ai[i]);
                    } else {
                        for (int k = i + 1; k < N; ++k) t -= ai[k] * x[k];
                        if (nounit) t /= ai[i];
                    }
                    x[i] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    if (conj_a) {
                        for (int k = 0; k < i; ++k) t -= conjq(ai[k]) * x[k];
                        if (nounit) t /= conjq(ai[i]);
                    } else {
                        for (int k = 0; k < i; ++k) t -= ai[k] * x[k];
                        if (nounit) t /= ai[i];
                    }
                    x[i] = t;
                }
            }
        }
    } else {
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        if (TR == 'N') {
            if (UPLO == 'L') {
                for (int i = 0; i < N; ++i) {
                    const int ix = kx + i * incx;
                    if (x[ix] != zero) {
                        if (nounit) x[ix] /= A_(i, i);
                        const T xi = x[ix];
                        for (int k = i + 1; k < N; ++k) x[kx + k * incx] -= xi * A_(k, i);
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    const int ix = kx + i * incx;
                    if (x[ix] != zero) {
                        if (nounit) x[ix] /= A_(i, i);
                        const T xi = x[ix];
                        for (int k = 0; k < i; ++k) x[kx + k * incx] -= xi * A_(k, i);
                    }
                }
            }
        } else {
            const int conj_a = (TR == 'C');
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[kx + i * incx];
                    for (int k = i + 1; k < N; ++k) {
                        const T aki = conj_a ? conjq(A_(k, i)) : A_(k, i);
                        t -= aki * x[kx + k * incx];
                    }
                    if (nounit) t /= (conj_a ? conjq(A_(i, i)) : A_(i, i));
                    x[kx + i * incx] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[kx + i * incx];
                    for (int k = 0; k < i; ++k) {
                        const T aki = conj_a ? conjq(A_(k, i)) : A_(k, i);
                        t -= aki * x[kx + k * incx];
                    }
                    if (nounit) t /= (conj_a ? conjq(A_(i, i)) : A_(i, i));
                    x[kx + i * incx] = t;
                }
            }
        }
    }
}

/* ── Block-parallel variant ─────────────────────────────────────────
 *
 * LAPACK-blocked algorithm: walk the diagonal in NB-wide blocks; for
 * each block, call the unblocked xtrsv_ on the small diagonal block,
 * then issue a parallel xgemv to apply that block's contribution to
 * the rest of x.
 *
 * At kind16 every scalar op is a libquadmath call (~100 ns), so the
 * unblocked reference is compute-bound and parallelism in the
 * trailing xgemv is what scales. Fork-join cost (~5 µs) is amortized
 * by the per-block trailing-update work (jb · mt quadmath multiplies,
 * each ~100 ns). The diagonal sub-solve stays sequential but is
 * bounded by NB² ops so the critical path is O(N·NB).
 *
 * Falls back to xtrsv_ for incx≠1 (parallelism not worth the
 * scattered-index xgemv) or N below 2·NB (overhead exceeds gain).
 */

extern void xgemv_(
    const char *trans,
    const int *m, const int *n,
    const T *alpha,
    const T *a, const int *lda,
    const T *x, const int *incx,
    const T *beta,
    T *y, const int *incy,
    size_t trans_len);

static int xtrsv_blocked_nb(void) {
    static int cached = 0;
    if (cached == 0) {
        const char *s = getenv("XTRSV_NB");
        int v = (s && *s) ? atoi(s) : 0;
        cached = (v > 0) ? v : XTRSV_BLOCKED_NB_DEFAULT;
    }
    return cached;
}

void xtrsv_blocked_(
    const char *uplo, const char *trans, const char *diag,
    const int *n_,
    const T *restrict a, const int *lda_,
    T *restrict x, const int *incx_,
    size_t uplo_len, size_t trans_len, size_t diag_len)
{
    const int N = *n_;
    const int lda = *lda_, incx = *incx_;
    const int nb = xtrsv_blocked_nb();
    const char UPLO = up(uplo);
    const char TR   = up(trans);

    if (N == 0) return;
    if (incx != 1 || N < 2 * nb) {
        xtrsv_(uplo, trans, diag, n_, a, lda_, x, incx_,
               uplo_len, trans_len, diag_len);
        return;
    }

    const T neg_one = -1.0Q + 0.0Qi;
    const T one_v   =  1.0Q + 0.0Qi;
    const char NN[1] = {'N'};
    const char TT[1] = {(TR == 'C') ? 'C' : 'T'};
    const int one_i = 1;

    if (TR == 'N') {
        if (UPLO == 'L') {
            /* Forward: solve A11 x1 = b1, then x2 -= A21 x1, repeat. */
            for (int j = 0; j < N; j += nb) {
                int jb = (N - j < nb) ? (N - j) : nb;
                xtrsv_(uplo, trans, diag, &jb, &A_(j, j), lda_,
                       &x[j], &one_i, uplo_len, trans_len, diag_len);
                int mt = N - j - jb;
                if (mt > 0) {
                    int j2 = j + jb;
                    xgemv_(NN, &mt, &jb, &neg_one, &A_(j2, j), lda_,
                           &x[j], &one_i, &one_v, &x[j2], &one_i, 1);
                }
            }
        } else {
            /* Backward: solve A22 x2 = b2, then x1 -= A12 x2, repeat. */
            int j = ((N - 1) / nb) * nb;
            while (j >= 0) {
                int jb = (N - j < nb) ? (N - j) : nb;
                xtrsv_(uplo, trans, diag, &jb, &A_(j, j), lda_,
                       &x[j], &one_i, uplo_len, trans_len, diag_len);
                if (j > 0) {
                    xgemv_(NN, &j, &jb, &neg_one, &A_(0, j), lda_,
                           &x[j], &one_i, &one_v, &x[0], &one_i, 1);
                }
                j -= nb;
            }
        }
    } else {
        /* TR='T' or 'C': op(A) x = b. For UPLO='L' we descend from the
         * last diagonal block (because op(A) is upper); for UPLO='U'
         * we ascend. Trailing update is xgemv with op = T/C on a
         * rectangular slice of A. */
        if (UPLO == 'L') {
            int j = ((N - 1) / nb) * nb;
            while (j >= 0) {
                int jb = (N - j < nb) ? (N - j) : nb;
                xtrsv_(uplo, trans, diag, &jb, &A_(j, j), lda_,
                       &x[j], &one_i, uplo_len, trans_len, diag_len);
                if (j > 0) {
                    /* x[0:j] -= op(A[j:j+jb, 0:j]) * x[j:j+jb].
                     * xgemv(op, M=jb, N=j) on submatrix &A_(j, 0). */
                    xgemv_(TT, &jb, &j, &neg_one, &A_(j, 0), lda_,
                           &x[j], &one_i, &one_v, &x[0], &one_i, 1);
                }
                j -= nb;
            }
        } else {
            for (int j = 0; j < N; j += nb) {
                int jb = (N - j < nb) ? (N - j) : nb;
                xtrsv_(uplo, trans, diag, &jb, &A_(j, j), lda_,
                       &x[j], &one_i, uplo_len, trans_len, diag_len);
                int mt = N - j - jb;
                if (mt > 0) {
                    int j2 = j + jb;
                    /* x[j2:N] -= op(A[j:j+jb, j2:N]) * x[j:j+jb] */
                    xgemv_(TT, &jb, &mt, &neg_one, &A_(j, j2), lda_,
                           &x[j], &one_i, &one_v, &x[j2], &one_i, 1);
                }
            }
        }
    }
}

#undef A_
