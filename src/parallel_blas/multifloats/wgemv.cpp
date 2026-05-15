/*
 * wgemv — multifloats complex DD general matrix-vector multiply.
 */

#include <cstddef>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {

#define WGEMV_OMP_MIN 64

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_cdd{ R{0.0, 0.0}, R{0.0, 0.0} };
const T one_cdd { R{1.0, 0.0}, R{0.0, 0.0} };

inline bool cdd_iszero(const T &x) {
    return x.re.limbs[0] == 0.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}
inline bool cdd_isone(const T &x) {
    return x.re.limbs[0] == 1.0 && x.re.limbs[1] == 0.0
        && x.im.limbs[0] == 0.0 && x.im.limbs[1] == 0.0;
}

inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }

} /* namespace */

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]

extern "C" void wgemv_(
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
    const char TR = up(trans);

    if (M == 0 || N == 0) return;

    const int leny = (TR == 'N') ? M : N;

    if (!cdd_isone(beta)) {
        int iy = (incy < 0) ? -(leny - 1) * incy : 0;
        for (int i = 0; i < leny; ++i) {
            if (cdd_iszero(beta)) y[iy] = zero_cdd;
            else                  y[iy] = cmul(y[iy], beta);
            iy += incy;
        }
    }
    if (cdd_iszero(alpha)) return;

    if (TR == 'N' && incx == 1 && incy == 1) {
#ifdef _OPENMP
        const int use_omp = (M >= WGEMV_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel if(use_omp)
        {
            int tid = 0, nt = 1;
            if (use_omp) { tid = omp_get_thread_num(); nt = omp_get_num_threads(); }
            const int i_lo = (static_cast<long long>(M) * tid) / nt;
            const int i_hi = (static_cast<long long>(M) * (tid + 1)) / nt;
            for (int j = 0; j < N; ++j) {
                const T xj = x[j];
                if (!cdd_iszero(xj)) {
                    const T t = cmul(alpha, xj);
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] = cadd(y[i], cmul(t, aj[i]));
                }
            }
        }
#else
        for (int j = 0; j < N; ++j) {
            const T xj = x[j];
            if (!cdd_iszero(xj)) {
                const T t = cmul(alpha, xj);
                const T *aj = &A_(0, j);
                for (int i = 0; i < M; ++i) y[i] = cadd(y[i], cmul(t, aj[i]));
            }
        }
#endif
    } else if ((TR == 'T' || TR == 'C') && incx == 1 && incy == 1) {
        const int conj_a = (TR == 'C');
#ifdef _OPENMP
        const int use_omp = (N >= WGEMV_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const T *aj = &A_(0, j);
            T s = zero_cdd;
            if (conj_a) {
                for (int i = 0; i < M; ++i) s = cadd(s, cmul(cconj(aj[i]), x[i]));
            } else {
                for (int i = 0; i < M; ++i) s = cadd(s, cmul(aj[i], x[i]));
            }
            y[j] = cadd(y[j], cmul(alpha, s));
        }
    } else {
        if (TR == 'N') {
            int jx = (incx < 0) ? -(N - 1) * incx : 0;
            for (int j = 0; j < N; ++j) {
                const T xj = x[jx];
                if (!cdd_iszero(xj)) {
                    const T t = cmul(alpha, xj);
                    int iy = (incy < 0) ? -(M - 1) * incy : 0;
                    for (int i = 0; i < M; ++i) {
                        y[iy] = cadd(y[iy], cmul(t, A_(i, j)));
                        iy += incy;
                    }
                }
                jx += incx;
            }
        } else {
            const int conj_a = (TR == 'C');
            int jy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int j = 0; j < N; ++j) {
                T s = zero_cdd;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                for (int i = 0; i < M; ++i) {
                    const T aij = conj_a ? cconj(A_(i, j)) : A_(i, j);
                    s = cadd(s, cmul(aij, x[ix]));
                    ix += incx;
                }
                y[jy] = cadd(y[jy], cmul(alpha, s));
                jy += incy;
            }
        }
    }
}

#undef A_
