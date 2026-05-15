/* msymv — multifloats real DD symmetric matrix-vector. */

#include <cstddef>
#include <cctype>
#include <multifloats.h>

namespace mf = multifloats;
using T = mf::float64x2;

namespace {
inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}
const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (T x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void msymv_(
    const char *uplo,
    const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *x, const int *incx_,
    const T *beta_,
    T *y, const int *incy_,
    std::size_t uplo_len)
{
    (void)uplo_len;
    const int N = *n_;
    const int lda = *lda_, incx = *incx_, incy = *incy_;
    const T alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);

    if (N == 0) return;

    if (!dd_isone(beta)) {
        int iy = (incy < 0) ? -(N - 1) * incy : 0;
        for (int i = 0; i < N; ++i) {
            if (dd_iszero(beta)) y[iy] = zero_dd;
            else                 y[iy] = y[iy] * beta;
            iy += incy;
        }
    }
    if (dd_iszero(alpha)) return;

    if (incx == 1 && incy == 1) {
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = zero_dd;
                const T *ai = &A_(0, i);
                y[i] = y[i] + temp1 * ai[i];
                for (int k = i + 1; k < N; ++k) {
                    y[k]  = y[k] + temp1 * ai[k];
                    temp2 = temp2 + ai[k] * x[k];
                }
                y[i] = y[i] + alpha * temp2;
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[i];
                T temp2 = zero_dd;
                const T *ai = &A_(0, i);
                for (int k = 0; k < i; ++k) {
                    y[k]  = y[k] + temp1 * ai[k];
                    temp2 = temp2 + ai[k] * x[k];
                }
                y[i] = y[i] + temp1 * ai[i] + alpha * temp2;
            }
        }
    } else {
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        int ky = (incy < 0) ? -(N - 1) * incy : 0;
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = zero_dd;
                y[ky + i * incy] = y[ky + i * incy] + temp1 * A_(i, i);
                for (int k = i + 1; k < N; ++k) {
                    y[ky + k * incy] = y[ky + k * incy] + temp1 * A_(k, i);
                    temp2 = temp2 + A_(k, i) * x[kx + k * incx];
                }
                y[ky + i * incy] = y[ky + i * incy] + alpha * temp2;
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = alpha * x[kx + i * incx];
                T temp2 = zero_dd;
                for (int k = 0; k < i; ++k) {
                    y[ky + k * incy] = y[ky + k * incy] + temp1 * A_(k, i);
                    temp2 = temp2 + A_(k, i) * x[kx + k * incx];
                }
                y[ky + i * incy] = y[ky + i * incy] + temp1 * A_(i, i) + alpha * temp2;
            }
        }
    }
}

#undef A_
