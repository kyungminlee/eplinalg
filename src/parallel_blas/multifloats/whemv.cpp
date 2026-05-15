/* whemv — multifloats Hermitian matrix-vector. */

#include <cstddef>
#include <cctype>
#include <multifloats.h>

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {
inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}
const R rzero{0.0, 0.0};
const T zero_cdd{ rzero, rzero };
const T one_cdd { R{1.0, 0.0}, rzero };

inline bool dd_iszero(R x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool cdd_iszero(const T &x) { return dd_iszero(x.re) && dd_iszero(x.im); }
inline bool cdd_isone (const T &x) {
    return x.re.limbs[0] == 1.0 && x.re.limbs[1] == 0.0 && dd_iszero(x.im);
}
inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void whemv_(
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

    if (!cdd_isone(beta)) {
        int iy = (incy < 0) ? -(N - 1) * incy : 0;
        for (int i = 0; i < N; ++i) {
            if (cdd_iszero(beta)) y[iy] = zero_cdd;
            else                  y[iy] = cmul(y[iy], beta);
            iy += incy;
        }
    }
    if (cdd_iszero(alpha)) return;

    if (incx == 1 && incy == 1) {
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = cmul(alpha, x[i]);
                T temp2 = zero_cdd;
                const T *ai = &A_(0, i);
                /* Diagonal: treat as real. */
                const T aii_re{ ai[i].re, rzero };
                y[i] = cadd(y[i], cmul(temp1, aii_re));
                for (int k = i + 1; k < N; ++k) {
                    y[k]  = cadd(y[k], cmul(temp1, ai[k]));
                    temp2 = cadd(temp2, cmul(cconj(ai[k]), x[k]));
                }
                y[i] = cadd(y[i], cmul(alpha, temp2));
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = cmul(alpha, x[i]);
                T temp2 = zero_cdd;
                const T *ai = &A_(0, i);
                for (int k = 0; k < i; ++k) {
                    y[k]  = cadd(y[k], cmul(temp1, ai[k]));
                    temp2 = cadd(temp2, cmul(cconj(ai[k]), x[k]));
                }
                const T aii_re{ ai[i].re, rzero };
                y[i] = cadd(y[i], cadd(cmul(temp1, aii_re), cmul(alpha, temp2)));
            }
        }
    } else {
        int kx = (incx < 0) ? -(N - 1) * incx : 0;
        int ky = (incy < 0) ? -(N - 1) * incy : 0;
        if (UPLO == 'L') {
            for (int i = 0; i < N; ++i) {
                const T temp1 = cmul(alpha, x[kx + i * incx]);
                T temp2 = zero_cdd;
                const T aii_re{ A_(i, i).re, rzero };
                y[ky + i * incy] = cadd(y[ky + i * incy], cmul(temp1, aii_re));
                for (int k = i + 1; k < N; ++k) {
                    y[ky + k * incy] = cadd(y[ky + k * incy], cmul(temp1, A_(k, i)));
                    temp2 = cadd(temp2, cmul(cconj(A_(k, i)), x[kx + k * incx]));
                }
                y[ky + i * incy] = cadd(y[ky + i * incy], cmul(alpha, temp2));
            }
        } else {
            for (int i = 0; i < N; ++i) {
                const T temp1 = cmul(alpha, x[kx + i * incx]);
                T temp2 = zero_cdd;
                for (int k = 0; k < i; ++k) {
                    y[ky + k * incy] = cadd(y[ky + k * incy], cmul(temp1, A_(k, i)));
                    temp2 = cadd(temp2, cmul(cconj(A_(k, i)), x[kx + k * incx]));
                }
                const T aii_re{ A_(i, i).re, rzero };
                y[ky + i * incy] = cadd(y[ky + i * incy], cadd(cmul(temp1, aii_re), cmul(alpha, temp2)));
            }
        }
    }
}

#undef A_
