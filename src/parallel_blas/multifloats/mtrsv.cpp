/* mtrsv — multifloats real DD triangular solve. */

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
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void mtrsv_(
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
    char TR = up(trans);
    if (TR == 'C') TR = 'T';
    const char DIAG = up(diag);
    const bool nounit = (DIAG != 'U');

    if (N == 0) return;

    if (incx == 1) {
        if (TR == 'N') {
            if (UPLO == 'L') {
                for (int i = 0; i < N; ++i) {
                    if (!dd_iszero(x[i])) {
                        if (nounit) x[i] = x[i] / A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = i + 1; k < N; ++k) x[k] = x[k] - xi * ai[k];
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    if (!dd_iszero(x[i])) {
                        if (nounit) x[i] = x[i] / A_(i, i);
                        const T xi = x[i];
                        const T *ai = &A_(0, i);
                        for (int k = 0; k < i; ++k) x[k] = x[k] - xi * ai[k];
                    }
                }
            }
        } else {
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    for (int k = i + 1; k < N; ++k) t = t - ai[k] * x[k];
                    if (nounit) t = t / ai[i];
                    x[i] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[i];
                    const T *ai = &A_(0, i);
                    for (int k = 0; k < i; ++k) t = t - ai[k] * x[k];
                    if (nounit) t = t / ai[i];
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
                    if (!dd_iszero(x[ix])) {
                        if (nounit) x[ix] = x[ix] / A_(i, i);
                        const T xi = x[ix];
                        for (int k = i + 1; k < N; ++k) x[kx + k * incx] = x[kx + k * incx] - xi * A_(k, i);
                    }
                }
            } else {
                for (int i = N - 1; i >= 0; --i) {
                    const int ix = kx + i * incx;
                    if (!dd_iszero(x[ix])) {
                        if (nounit) x[ix] = x[ix] / A_(i, i);
                        const T xi = x[ix];
                        for (int k = 0; k < i; ++k) x[kx + k * incx] = x[kx + k * incx] - xi * A_(k, i);
                    }
                }
            }
        } else {
            if (UPLO == 'L') {
                for (int i = N - 1; i >= 0; --i) {
                    T t = x[kx + i * incx];
                    for (int k = i + 1; k < N; ++k) t = t - A_(k, i) * x[kx + k * incx];
                    if (nounit) t = t / A_(i, i);
                    x[kx + i * incx] = t;
                }
            } else {
                for (int i = 0; i < N; ++i) {
                    T t = x[kx + i * incx];
                    for (int k = 0; k < i; ++k) t = t - A_(k, i) * x[kx + k * incx];
                    if (nounit) t = t / A_(i, i);
                    x[kx + i * incx] = t;
                }
            }
        }
    }
}

#undef A_
