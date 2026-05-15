/*
 * mgemv — multifloats real DD general matrix-vector multiply.
 */

#include <cstddef>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using T = mf::float64x2;

namespace {

#define MGEMV_OMP_MIN 64

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (T x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }

} /* namespace */

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void mgemv_(
    const char *trans,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *x, const int *incx_,
    const T *beta_,
    T *y, const int *incy_,
    std::size_t trans_len)
{
    (void)trans_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, incx = *incx_, incy = *incy_;
    const T alpha = *alpha_, beta = *beta_;
    char TR = up(trans);
    if (TR == 'C') TR = 'T';

    if (M == 0 || N == 0) return;

    const int leny = (TR == 'N') ? M : N;

    if (!dd_isone(beta)) {
        int iy = (incy < 0) ? -(leny - 1) * incy : 0;
        for (int i = 0; i < leny; ++i) {
            if (dd_iszero(beta)) y[iy] = zero_dd;
            else                 y[iy] = y[iy] * beta;
            iy += incy;
        }
    }
    if (dd_iszero(alpha)) return;

    if (TR == 'N' && incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (M >= MGEMV_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel if(use_omp)
        {
            int tid = 0, nt = 1;
            if (use_omp) { tid = omp_get_thread_num(); nt = omp_get_num_threads(); }
            const int i_lo = (static_cast<long long>(M) * tid) / nt;
            const int i_hi = (static_cast<long long>(M) * (tid + 1)) / nt;
            for (int j = 0; j < N; ++j) {
                const T xj = x[j];
                if (!dd_iszero(xj)) {
                    const T t = alpha * xj;
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] = y[i] + t * aj[i];
                }
            }
        }
#else
        for (int j = 0; j < N; ++j) {
            const T xj = x[j];
            if (!dd_iszero(xj)) {
                const T t = alpha * xj;
                const T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) y[i] = y[i] + t * aj[i];
            }
        }
#endif
    } else if (TR != 'N' && incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (N >= MGEMV_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const T *aj = &A_(0, j);
            T s = zero_dd;
            for (int i = 0; i < M; ++i) s = s + aj[i] * x[i];
            y[j] = y[j] + alpha * s;
        }
    } else {
        if (TR == 'N') {
            int jx = (incx < 0) ? -(N - 1) * incx : 0;
            for (int j = 0; j < N; ++j) {
                const T xj = x[jx];
                if (!dd_iszero(xj)) {
                    const T t = alpha * xj;
                    int iy = (incy < 0) ? -(M - 1) * incy : 0;
                    for (int i = 0; i < M; ++i) {
                        y[iy] = y[iy] + t * A_(i, j);
                        iy += incy;
                    }
                }
                jx += incx;
            }
        } else {
            int jy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int j = 0; j < N; ++j) {
                T s = zero_dd;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                for (int i = 0; i < M; ++i) {
                    s = s + A_(i, j) * x[ix];
                    ix += incx;
                }
                y[jy] = y[jy] + alpha * s;
                jy += incy;
            }
        }
    }
}

#undef A_
