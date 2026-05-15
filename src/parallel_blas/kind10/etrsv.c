/*
 * etrsv — kind10 (REAL(KIND=10)) triangular solve.
 *   A x = b           (TRANS='N')
 *   Aᵀ x = b          (TRANS='T'/'C')
 * where A is N×N triangular (UPLO, DIAG). x overwrites b in-place.
 *
 * Inherently serial in i (each x[i] depends on earlier-solved x[k]),
 * so no OMP. Mirrors the Netlib reference with restrict and stride-1
 * column access of A where possible.
 */

#include <stddef.h>
#include <ctype.h>

typedef long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void etrsv_(
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
    char TR = up(trans);
    if (TR == 'C') TR = 'T';
    const char DIAG = up(diag);
    const int nounit = (DIAG != 'U');

    if (N == 0) return;

    const T zero = 0.0L;

    if (incx == 1) {
        if (TR == 'N') {
            if (UPLO == 'L') {
                /* Forward substitution: x[i] = (b[i] - sum_{k<i} A(i,k) x[k]) / A(i,i). */
                for (int i = 0; i < N; ++i) {
                    if (x[i] != zero) {
                        if (nounit) x[i] /= A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        /* AXPY-style update of x[i+1..N-1] using column i of A. */
                        for (int k = i + 1; k < N; ++k) x[k] -= xi * ai[k];
                    }
                }
            } else {
                /* UPLO='U': back substitution iterates i backward.
                 * x[i] = (b[i] - sum_{k>i} A(i,k) x[k]) / A(i,i). */
                for (int i = N - 1; i >= 0; --i) {
                    if (x[i] != zero) {
                        if (nounit) x[i] /= A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = 0; k < i; ++k) x[k] -= xi * ai[k];
                    }
                }
            }
        } else {  /* TRANS = 'T': solve Aᵀ x = b. */
            if (UPLO == 'L') {
                /* Lower-stored A; Aᵀ is upper. Iterate i backward.
                 * x[i] = (b[i] - sum_{k>i} A(k,i) x[k]) / A(i,i). */
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    for (int k = i + 1; k < N; ++k) t -= ai[k] * x[k];
                    if (nounit) t /= ai[i];
                    x[i] = t;
                }
            } else {
                /* UPLO='U': iterate i forward.
                 * x[i] = (b[i] - sum_{k<i} A(k,i) x[k]) / A(i,i). */
                for (int i = 0; i < N; ++i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    for (int k = 0; k < i; ++k) t -= ai[k] * x[k];
                    if (nounit) t /= ai[i];
                    x[i] = t;
                }
            }
        }
    } else {
        /* General-stride fallback. */
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
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[kx + i * incx];
                    for (int k = i + 1; k < N; ++k) t -= A_(k, i) * x[kx + k * incx];
                    if (nounit) t /= A_(i, i);
                    x[kx + i * incx] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[kx + i * incx];
                    for (int k = 0; k < i; ++k) t -= A_(k, i) * x[kx + k * incx];
                    if (nounit) t /= A_(i, i);
                    x[kx + i * incx] = t;
                }
            }
        }
    }
}

#undef A_
