/* mger — multifloats real DD rank-1 update. */

#include <cstddef>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using T = mf::float64x2;

namespace {
#define MGER_OMP_MIN 64
const T zero_dd{0.0, 0.0};
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
}

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void mger_(
    const int *m_, const int *n_,
    const T *alpha_,
    const T *x, const int *incx_,
    const T *y, const int *incy_,
    T *a, const int *lda_)
{
    const int M = *m_, N = *n_;
    const int incx = *incx_, incy = *incy_, lda = *lda_;
    const T alpha = *alpha_;

    if (M == 0 || N == 0 || dd_iszero(alpha)) return;

    if (incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (N >= MGER_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const T yj = y[j];
            if (!dd_iszero(yj)) {
                const T t = alpha * yj;
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) aj[i] = aj[i] + t * x[i];
            }
        }
    } else {
        int jy = (incy < 0) ? -(N - 1) * incy : 0;
        for (int j = 0; j < N; ++j) {
            const T yj = y[jy];
            if (!dd_iszero(yj)) {
                const T t = alpha * yj;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) {
                    aj[i] = aj[i] + t * x[ix];
                    ix += incx;
                }
            }
            jy += incy;
        }
    }
}

#undef A_
