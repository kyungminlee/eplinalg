/*
 * ygemv — kind10 complex (`_Complex long double`) general matrix-vector.
 *
 *   y := alpha · A · x + beta · y            (TRANS='N')   A is M×N
 *   y := alpha · Aᵀ · x + beta · y           (TRANS='T')   A is M×N, y is N
 *   y := alpha · Aᴴ · x + beta · y           (TRANS='C')   A is M×N, y is N
 *
 * Reference ZGEMV with J-axis unroll by 2 on the N path (halves y
 * memory traffic) and K-axis unroll by 2 on the T/C dot-product paths
 * (two independent fadd chains hide x87 latency on the cmul-accumulate
 * pattern).
 *
 * Partitioning strategy: thread owns a slice of the output vector y.
 * For TRANS='N' that's rows of A; for TRANS='T'/'C' that's columns of A.
 * All A accesses inside a thread's slice are column-major stride-1.
 */

#include <stddef.h>
#include <ctype.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define YGEMV_OMP_MIN 64

typedef _Complex long double T;

static inline char up(const char *p) {
    return (char)toupper((unsigned char)*p);
}

static inline T cconj(T z) { return ~z; }

static const T ZERO = 0.0L + 0.0Li;
static const T ONE  = 1.0L + 0.0Li;

#define A_(i, j)  a[(size_t)(j) * lda + (i)]

void ygemv_(
    const char *trans,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *restrict a, const int *lda_,
    const T *restrict x, const int *incx_,
    const T *beta_,
    T *restrict y, const int *incy_,
    size_t trans_len)
{
    (void)trans_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, incx = *incx_, incy = *incy_;
    const T alpha = *alpha_, beta = *beta_;
    const char TR = up(trans);

    if (M == 0 || N == 0) return;

    const int leny = (TR == 'N') ? M : N;

    /* β-scale y. */
    if (beta != ONE) {
        if (incy == 1) {
            if (beta == ZERO) for (int i = 0; i < leny; ++i) y[i] = ZERO;
            else              for (int i = 0; i < leny; ++i) y[i] *= beta;
        } else {
            int iy = (incy < 0) ? -(leny - 1) * incy : 0;
            for (int i = 0; i < leny; ++i) {
                if (beta == ZERO) y[iy] = ZERO;
                else              y[iy] *= beta;
                iy += incy;
            }
        }
    }

    if (alpha == ZERO) return;

    if (TR == 'N') {
        if (incx == 1 && incy == 1) {
#ifdef _OPENMP
            const int use_omp = (M >= YGEMV_OMP_MIN && omp_get_max_threads() > 1);
#endif
            int i_lo = 0, i_hi = M;
            (void)i_lo; (void)i_hi;
#ifdef _OPENMP
            #pragma omp parallel if(use_omp) firstprivate(i_lo, i_hi)
            {
                if (use_omp) {
                    const int tid = omp_get_thread_num();
                    const int nt  = omp_get_num_threads();
                    i_lo = ((long long)M * tid) / nt;
                    i_hi = ((long long)M * (tid + 1)) / nt;
                }
#endif
                int j = 0;
                for (; j + 1 < N; j += 2) {
                    const T t0 = alpha * x[j];
                    const T t1 = alpha * x[j + 1];
                    const T *a0 = &A_(0, j);
                    const T *a1 = &A_(0, j + 1);
                    for (int i = i_lo; i < i_hi; ++i) {
                        y[i] += t0 * a0[i] + t1 * a1[i];
                    }
                }
                for (; j < N; ++j) {
                    const T t = alpha * x[j];
                    const T *aj = &A_(0, j);
                    for (int i = i_lo; i < i_hi; ++i) y[i] += t * aj[i];
                }
#ifdef _OPENMP
            }
#endif
        } else {
            int jx = (incx < 0) ? -(N - 1) * incx : 0;
            for (int j = 0; j < N; ++j) {
                const T xj = x[jx];
                if (xj != ZERO) {
                    const T t = alpha * xj;
                    int iy = (incy < 0) ? -(M - 1) * incy : 0;
                    for (int i = 0; i < M; ++i) {
                        y[iy] += t * A_(i, j);
                        iy += incy;
                    }
                }
                jx += incx;
            }
        }
    } else {
        /* TRANS='T' or 'C': y[j] += α · Σ_i A(i,j)[, conj] · x(i).
         * Single-acc dot product (complex inner-loop is x87-stack-heavy;
         * K-unroll with split accumulators regressed on similar paths
         * in ygemm — keep single accumulator). */
        const int conj_a = (TR == 'C');
        if (incx == 1 && incy == 1) {
#ifdef _OPENMP
            const int use_omp = (N >= YGEMV_OMP_MIN && omp_get_max_threads() > 1);
            #pragma omp parallel for if(use_omp) schedule(static)
#endif
            for (int j = 0; j < N; ++j) {
                const T *aj = &A_(0, j);
                T s = ZERO;
                if (conj_a) {
                    for (int i = 0; i < M; ++i) s += cconj(aj[i]) * x[i];
                } else {
                    for (int i = 0; i < M; ++i) s += aj[i] * x[i];
                }
                y[j] += alpha * s;
            }
        } else {
            int jy = (incy < 0) ? -(N - 1) * incy : 0;
            for (int j = 0; j < N; ++j) {
                T s = ZERO;
                int ix = (incx < 0) ? -(M - 1) * incx : 0;
                for (int i = 0; i < M; ++i) {
                    s += (conj_a ? cconj(A_(i, j)) : A_(i, j)) * x[ix];
                    ix += incx;
                }
                y[jy] += alpha * s;
                jy += incy;
            }
        }
    }
}

#undef A_
