/*
 * yhemv — kind10 complex Hermitian matrix-vector multiply.
 *   y := alpha · A · x + beta · y    where A is N×N Hermitian
 *
 * Same two-pass pattern as esymv. Hermitian twist: the cross
 * reflection conjugates A. Diagonal of A is treated as real (matches
 * Netlib ZHEMV — uses DBLE(A(I,I))).
 */

#include <stddef.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define YHEMV_OMP_MIN 64

typedef _Complex long double T;
static const T ZERO = 0.0L + 0.0Li;
static const T ONE  = 1.0L + 0.0Li;
static inline T cconj(T z) { return ~z; }

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void yhemv_(
    const char *uplo,
    const int *n_,
    const T *alpha_,
    const T *restrict a, const int *lda_,
    const T *restrict x, const int *incx_,
    const T *beta_,
    T *restrict y, const int *incy_,
    size_t uplo_len)
{
    (void)uplo_len;
    const int N = *n_;
    const int lda = *lda_, incx = *incx_, incy = *incy_;
    const T alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);

    if (N == 0) return;

    if (beta != ONE) {
        if (incy == 1) {
            if (beta == ZERO) for (int i = 0; i < N; ++i) y[i] = ZERO;
            else              for (int i = 0; i < N; ++i) y[i] *= beta;
        } else {
            int iy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int i = 0; i < N; ++i) {
                if (beta == ZERO) y[iy] = ZERO;
                else              y[iy] *= beta;
                iy += incy;
            }
        }
    }

    if (alpha == ZERO) return;

    if (incx == 1 && incy == 1) {
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = ZERO;
                const T *ai = &A_(0, i);
                y[i] += temp1 * __real__ ai[i];
                for (int k = i + 1; k < N; ++k) {
                    y[k]  += temp1 * ai[k];
                    temp2 += cconj(ai[k]) * x[k];
                }
                y[i] += alpha * temp2;
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = ZERO;
                const T *ai = &A_(0, i);
                for (int k = 0; k < i; ++k) {
                    y[k]  += temp1 * ai[k];
                    temp2 += cconj(ai[k]) * x[k];
                }
                y[i] += temp1 * __real__ ai[i] + alpha * temp2;
            }
        }
    } else {
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        int ky = (incy < 0) ? -(N - 1) * incy : 0;
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = ZERO;
                y[ky + i * incy] += temp1 * __real__ A_(i, i);
                for (int k = i + 1; k < N; ++k) {
                    y[ky + k * incy] += temp1 * A_(k, i);
                    temp2 += cconj(A_(k, i)) * x[kx + k * incx];
                }
                y[ky + i * incy] += alpha * temp2;
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = ZERO;
                for (int k = 0; k < i; ++k) {
                    y[ky + k * incy] += temp1 * A_(k, i);
                    temp2 += cconj(A_(k, i)) * x[kx + k * incx];
                }
                y[ky + i * incy] += temp1 * __real__ A_(i, i) + alpha * temp2;
            }
        }
    }
}

#undef A_
