/*
 * esymv — kind10 (REAL(KIND=10)) symmetric matrix-vector multiply.
 *   y := alpha · A · x + beta · y    where A is N×N symmetric
 *
 * Uses Netlib DSYMV's two-pass pattern: for each i,
 *   temp1 = alpha · x(i)   (contributes to y(k) for k!=i via A column reads)
 *   temp2 = sum_k A(k,i) · x(k)   (dot-product accumulator)
 *   y(i) += temp1 · A(i,i) + alpha · temp2
 * Stride-1 walks of A by columns; same direction-flip trick as the
 * symm/hemm diagonal kernel.
 */

#include <stddef.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define ESYMV_OMP_MIN 64

typedef long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void esymv_(
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

    const T zero = 0.0L, one = 1.0L;

    if (beta != one) {
        if (incy == 1) {
            if (beta == zero) for (int i = 0; i < N; ++i) y[i] = zero;
            else              for (int i = 0; i < N; ++i) y[i] *= beta;
        } else {
            int iy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int i = 0; i < N; ++i) {
                if (beta == zero) y[iy] = zero;
                else              y[iy] *= beta;
                iy += incy;
            }
        }
    }

    if (alpha == zero) return;

    /* The unit-stride path: stride-1 column walks of A. */
    if (incx == 1 && incy == 1) {
        if (UPLO == 'L') {
            /* Iterate i forward; the inner k loop covers k = i..N-1
             * (stored lower triangle). Uses A_(k, i) (stride-1 in k). */
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = zero;
                const T *ai = &A_(0, i);
                y[i] += temp1 * ai[i];
                for (int k = i + 1; k < N; ++k) {
                    y[k]  += temp1 * ai[k];
                    temp2 += ai[k] * x[k];
                }
                y[i] += alpha * temp2;
            }
        } else {
            /* UPLO='U': iterate i forward; inner k = 0..i-1
             * (stored upper triangle). */
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = zero;
                const T *ai = &A_(0, i);
                for (int k = 0; k < i; ++k) {
                    y[k]  += temp1 * ai[k];
                    temp2 += ai[k] * x[k];
                }
                y[i] += temp1 * ai[i] + alpha * temp2;
            }
        }
    } else {
        /* General-stride fallback: faithful to Netlib reference. */
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        int ky = (incy < 0) ? -(N - 1) * incy : 0;
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = zero;
                y[ky + i * incy] += temp1 * A_(i, i);
                for (int k = i + 1; k < N; ++k) {
                    y[ky + k * incy] += temp1 * A_(k, i);
                    temp2 += A_(k, i) * x[kx + k * incx];
                }
                y[ky + i * incy] += alpha * temp2;
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = zero;
                for (int k = 0; k < i; ++k) {
                    y[ky + k * incy] += temp1 * A_(k, i);
                    temp2 += A_(k, i) * x[kx + k * incx];
                }
                y[ky + i * incy] += temp1 * A_(i, i) + alpha * temp2;
            }
        }
    }
}

#undef A_
