/* wgerc — multifloats complex DD conjugated rank-1. */

#include <cstddef>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {
#define WGERC_OMP_MIN 64
const T zero_cdd{ R{0.0, 0.0}, R{0.0, 0.0} };
inline bool cdd_iszero(const T &x) {
    return x.re.limbs[0] == 0.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}
inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void wgerc_(
    const int *m_, const int *n_,
    const T *alpha_,
    const T *x, const int *incx_,
    const T *y, const int *incy_,
    T *a, const int *lda_)
{
    const int M = *m_, N = *n_;
    const int incx = *incx_, incy = *incy_, lda = *lda_;
    const T alpha = *alpha_;

    if (M == 0 || N == 0 || cdd_iszero(alpha)) return;

    if (incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (N >= WGERC_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const T yj = cconj(y[j]);
            if (!cdd_iszero(yj)) {
                const T t = cmul(alpha, yj);
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) aj[i] = cadd(aj[i], cmul(t, x[i]));
            }
        }
    } else {
        int jy = (incy < 0) ? -(N - 1) * incy : 0;
        for (int j = 0; j < N; ++j) {
            const T yj = cconj(y[jy]);
            if (!cdd_iszero(yj)) {
                const T t = cmul(alpha, yj);
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) {
                    aj[i] = cadd(aj[i], cmul(t, x[ix]));
                    ix += incx;
                }
            }
            jy += incy;
        }
    }
}

#undef A_
