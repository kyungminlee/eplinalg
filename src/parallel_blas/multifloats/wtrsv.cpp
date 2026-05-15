/* wtrsv — multifloats complex DD triangular solve. */

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
const T zero_cdd{ R{0.0, 0.0}, R{0.0, 0.0} };
inline bool cdd_iszero(const T &x) {
    return x.re.limbs[0] == 0.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}
inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T csub(T const &a, T const &b) { return T{ a.re - b.re, a.im - b.im }; }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }
/* Complex DD division via the standard formula. */
inline T cdiv(T const &a, T const &b) {
    const R d = b.re * b.re + b.im * b.im;
    const R inv_d = R{1.0, 0.0} / d;
    return T{ (a.re * b.re + a.im * b.im) * inv_d,
              (a.im * b.re - a.re * b.im) * inv_d };
}
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void wtrsv_(
    const char *uplo, const char *trans, const char *diag,
    const int *n_,
    const T *a, const int *lda_,
    T *x, const int *incx_,
    std::size_t uplo_len, std::size_t trans_len, std::size_t diag_len)
{
    (void)uplo_len; (void)trans_len; (void)diag_len;
    const int N = *n_;
    const int lda = *lda_, incx = *incx_;
    const char UPLO = up(uplo);
    const char TR   = up(trans);
    const char DIAG = up(diag);
    const bool nounit = (DIAG != 'U');

    if (N == 0) return;

    if (incx == 1) {
        if (TR == 'N') {
            if (UPLO == 'L') {
                for (int i = 0; i < N; ++i) {
                    if (!cdd_iszero(x[i])) {
                        if (nounit) x[i] = cdiv(x[i], A_(i, i));
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = i + 1; k < N; ++k) x[k] = csub(x[k], cmul(xi, ai[k]));
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    if (!cdd_iszero(x[i])) {
                        if (nounit) x[i] = cdiv(x[i], A_(i, i));
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = 0; k < i; ++k) x[k] = csub(x[k], cmul(xi, ai[k]));
                    }
                }
            }
        } else {
            const bool conj_a = (TR == 'C');
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    if (conj_a) {
                        for (int k = i + 1; k < N; ++k) t = csub(t, cmul(cconj(ai[k]), x[k]));
                        if (nounit) t = cdiv(t, cconj(ai[i]));
                    } else {
                        for (int k = i + 1; k < N; ++k) t = csub(t, cmul(ai[k], x[k]));
                        if (nounit) t = cdiv(t, ai[i]);
                    }
                    x[i] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    if (conj_a) {
                        for (int k = 0; k < i; ++k) t = csub(t, cmul(cconj(ai[k]), x[k]));
                        if (nounit) t = cdiv(t, cconj(ai[i]));
                    } else {
                        for (int k = 0; k < i; ++k) t = csub(t, cmul(ai[k], x[k]));
                        if (nounit) t = cdiv(t, ai[i]);
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
                    if (!cdd_iszero(x[ix])) {
                        if (nounit) x[ix] = cdiv(x[ix], A_(i, i));
                        const T xi = x[ix];
                        for (int k = i + 1; k < N; ++k) x[kx + k * incx] = csub(x[kx + k * incx], cmul(xi, A_(k, i)));
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    const int ix = kx + i * incx;
                    if (!cdd_iszero(x[ix])) {
                        if (nounit) x[ix] = cdiv(x[ix], A_(i, i));
                        const T xi = x[ix];
                        for (int k = 0; k < i; ++k) x[kx + k * incx] = csub(x[kx + k * incx], cmul(xi, A_(k, i)));
                    }
                }
            }
        } else {
            const bool conj_a = (TR == 'C');
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[kx + i * incx];
                    for (int k = i + 1; k < N; ++k) {
                        const T aki = conj_a ? cconj(A_(k, i)) : A_(k, i);
                        t = csub(t, cmul(aki, x[kx + k * incx]));
                    }
                    if (nounit) t = cdiv(t, (conj_a ? cconj(A_(i, i)) : A_(i, i)));
                    x[kx + i * incx] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[kx + i * incx];
                    for (int k = 0; k < i; ++k) {
                        const T aki = conj_a ? cconj(A_(k, i)) : A_(k, i);
                        t = csub(t, cmul(aki, x[kx + k * incx]));
                    }
                    if (nounit) t = cdiv(t, (conj_a ? cconj(A_(i, i)) : A_(i, i)));
                    x[kx + i * incx] = t;
                }
            }
        }
    }
}

#undef A_
