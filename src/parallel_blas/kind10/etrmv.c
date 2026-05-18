/*
 * etrmv — kind10 (REAL(KIND=10)) triangular matrix-vector.
 *   x := A · x         (TRANS='N')
 *   x := Aᵀ · x        (TRANS='T'/'C')
 * A is N×N triangular (UPLO, DIAG). x updated in-place.
 *
 * Inherently serial over j (each j writes to x[j] or reads x[i] from
 * a region that earlier j's modified). No OMP. Netlib reference +
 * restrict + stride-1 column access.
 */

#include <stddef.h>
#include <ctype.h>

typedef long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void etrmv_(
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
                /* j backward: x[i] for i>j updated by temp=x[j]; then scale x[j].
                 * Inner walks backward (i = N-1..j+1) to match Fortran
                 * etrmv.f (DO 50 I = N,J+1,-1). Sub-class C / Rule 21. */
                for (int j = N - 1; j >= 0; --j) {
                    const T temp = x[j];
                    if (temp != zero) {
                        const T *aj = &A_(0, j);
                        for (int i = N - 1; i > j; --i) x[i] += temp * aj[i];
                    }
                    if (nounit) x[j] *= A_(j, j);
                }
            } else {
                /* UPLO='U', j forward: x[i] for i<j updated by temp; then scale. */
                for (int j = 0; j < N; ++j) {
                    const T temp = x[j];
                    if (temp != zero) {
                        const T *aj = &A_(0, j);
                        for (int i = 0; i < j; ++i) x[i] += temp * aj[i];
                    }
                    if (nounit) x[j] *= A_(j, j);
                }
            }
        } else {  /* TRANS = 'T' */
            if (UPLO == 'L') {
                /* j forward: dot product over i>j into x[j]. */
                for (int j = 0; j < N; ++j) {
                    T temp = x[j];
                    if (nounit) temp *= A_(j, j);
                    const T *aj = &A_(0, j);
                    /* 2-chain dot product (x87 latency-hiding). */
                    T s0 = zero, s1 = zero;
                    int i = j + 1;
                    for (; i + 1 < N; i += 2) {
                        s0 += aj[i]     * x[i];
                        s1 += aj[i + 1] * x[i + 1];
                    }
                    T s = s0 + s1;
                    for (; i < N; ++i) s += aj[i] * x[i];
                    x[j] = temp + s;
                }
            } else {
                /* UPLO='U', j backward: dot over i<j into x[j].
                 * Inner walks backward (i = j-1..0) to match the
                 * Fortran reference (DO 90 I = J-1,1,-1). 2-chain
                 * unroll preserved: descend in pairs. Sub-class D /
                 * Rule 21 — even though the current forward-2-chain
                 * already beats migrated at measured N because of
                 * x87 latency hiding, the backward walk keeps the
                 * direction consistent with the Fortran reference. */
                for (int j = N - 1; j >= 0; --j) {
                    T temp = x[j];
                    if (nounit) temp *= A_(j, j);
                    const T *aj = &A_(0, j);
                    T s0 = zero, s1 = zero;
                    int i = j - 1;
                    for (; i - 1 >= 0; i -= 2) {
                        s0 += aj[i]     * x[i];
                        s1 += aj[i - 1] * x[i - 1];
                    }
                    T s = s0 + s1;
                    for (; i >= 0; --i) s += aj[i] * x[i];
                    x[j] = temp + s;
                }
            }
        }
    } else {
        /* General-stride fallback. */
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        if (TR == 'N') {
            if (UPLO == 'L') {
                /* Inner walks backward to match Fortran etrmv.f
                 * (DO 70 I = N,J+1,-1). Sub-class C / Rule 21. */
                for (int j = N - 1; j >= 0; --j) {
                    const T temp = x[kx + j * incx];
                    if (temp != zero) {
                        for (int i = N - 1; i > j; --i) x[kx + i * incx] += temp * A_(i, j);
                    }
                    if (nounit) x[kx + j * incx] *= A_(j, j);
                }
            } else {
                for (int j = 0; j < N; ++j) {
                    const T temp = x[kx + j * incx];
                    if (temp != zero) {
                        for (int i = 0; i < j; ++i) x[kx + i * incx] += temp * A_(i, j);
                    }
                    if (nounit) x[kx + j * incx] *= A_(j, j);
                }
            }
        } else {
            if (UPLO == 'L') {
                for (int j = 0; j < N; ++j) {
                    T temp = x[kx + j * incx];
                    if (nounit) temp *= A_(j, j);
                    for (int i = j + 1; i < N; ++i) temp += A_(i, j) * x[kx + i * incx];
                    x[kx + j * incx] = temp;
                }
            } else {
                /* Inner walks backward to match Fortran reference
                 * (DO 110 I = J-1,1,-1). Sub-class D / Rule 21. */
                for (int j = N - 1; j >= 0; --j) {
                    T temp = x[kx + j * incx];
                    if (nounit) temp *= A_(j, j);
                    for (int i = j - 1; i >= 0; --i) temp += A_(i, j) * x[kx + i * incx];
                    x[kx + j * incx] = temp;
                }
            }
        }
    }
}

#undef A_
