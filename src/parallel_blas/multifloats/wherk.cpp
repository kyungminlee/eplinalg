/*
 * wherk — multifloats complex (DD) Hermitian rank-k.
 * alpha/beta REAL, A/C complex. Diagonal of C stays real.
 * Blocked: scalar diagonal (with real-diagonal accumulation) +
 * wgemm trailing with conjugate transpose.
 */

#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using R = mf::float64x2;
using T = mf::complex64x2;

namespace {

#define WHERK_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int herk_nb(void) { if (g_nb == 0) g_nb = env_int("WHERK_NB", 64); return g_nb; }

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const R rzero{0.0, 0.0};
const R rone {1.0, 0.0};
const T czero{ rzero, rzero };
const T cone { rone, rzero };

inline bool dd_iszero(R x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (R x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }
inline bool cdd_iszero(const T &x) { return dd_iszero(x.re) && dd_iszero(x.im); }

inline T cmul(T const &a, T const &b) {
    return T{ a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re };
}
inline T cadd(T const &a, T const &b) { return T{ a.re + b.re, a.im + b.im }; }
inline T cconj(T const &a) { return T{ a.re, R{-a.im.limbs[0], -a.im.limbs[1]} }; }
inline T rcmul(R const &r, T const &z) { return T{ r * z.re, r * z.im }; }

extern "C" void wgemm_(
    const char *transa, const char *transb,
    const int *m, const int *n, const int *k,
    const T *alpha,
    const T *a, const int *lda,
    const T *b, const int *ldb,
    const T *beta,
    T *c, const int *ldc,
    std::size_t transa_len, std::size_t transb_len);

#define A_(i, j)  a[static_cast<std::size_t>(j) * lda + (i)]
#define C_(i, j)  c[static_cast<std::size_t>(j) * ldc + (i)]

void herk_diag_add(int jc, int jb, int K, R alpha,
                   const T *a, int lda, T *c, int ldc,
                   char UPLO, char TR_c)
{
    if (TR_c == 'N') {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            for (int l = 0; l < K; ++l) {
                const T ajl = A_(j, l);
                if (!cdd_iszero(ajl)) {
                    const T t = rcmul(alpha, cconj(ajl));
                    for (int i = i_lo; i < i_hi; ++i) {
                        T prod = cmul(t, A_(i, l));
                        if (i == j) cj[i] = T{ cj[i].re + prod.re, cj[i].im };
                        else        cj[i] = cadd(cj[i], prod);
                    }
                }
            }
        }
    } else {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            const T *Aj = a + static_cast<std::size_t>(j) * lda;
            for (int i = i_lo; i < i_hi; ++i) {
                const T *Ai = a + static_cast<std::size_t>(i) * lda;
                T s = czero;
                for (int l = 0; l < K; ++l) s = cadd(s, cmul(cconj(Ai[l]), Aj[l]));
                T as = rcmul(alpha, s);
                if (i == j) cj[i] = T{ cj[i].re + as.re, cj[i].im };
                else        cj[i] = cadd(cj[i], as);
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void wherk_(
    const char *uplo, const char *trans,
    const int *n_, const int *k_,
    const R *alpha_,
    const T *a, const int *lda_,
    const R *beta_,
    T *c, const int *ldc_,
    std::size_t uplo_len, std::size_t trans_len)
{
    (void)uplo_len; (void)trans_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldc = *ldc_;
    const R alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);
    const char TR_c = up(trans);

    if (N == 0) return;

    const T alpha_c = T{ alpha, rzero };
    const char NN[1] = {'N'};
    const char CN[1] = {'C'};

    if (dd_iszero(alpha) || K == 0) {
        if (dd_isone(beta)) {
            for (int j = 0; j < N; ++j) {
                T *cjj = &c[static_cast<std::size_t>(j) * ldc + j];
                cjj->im = rzero;
            }
            return;
        }
#ifdef _OPENMP
        const bool use_omp = (N >= WHERK_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (dd_iszero(beta)) for (int i = i_lo; i < i_hi; ++i) cj[i] = czero;
            else {
                for (int i = i_lo; i < i_hi; ++i) {
                    if (i == j) cj[i] = T{ beta * cj[i].re, rzero };
                    else        cj[i] = rcmul(beta, cj[i]);
                }
            }
        }
        return;
    }

    const int nb = herk_nb();

#ifdef _OPENMP
    const bool use_omp = (N >= WHERK_OMP_MIN && omp_get_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
    for (int jc = 0; jc < N; jc += nb) {
        const int jb = (N - jc < nb) ? (N - jc) : nb;

        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (dd_iszero(beta)) {
                for (int i = i_lo; i < i_hi; ++i) cj[i] = czero;
            } else if (!dd_isone(beta)) {
                for (int i = i_lo; i < i_hi; ++i) {
                    if (i == j) cj[i] = T{ beta * cj[i].re, rzero };
                    else        cj[i] = rcmul(beta, cj[i]);
                }
            } else {
                cj[j].im = rzero;
            }
        }

        herk_diag_add(jc, jb, K, alpha, a, lda, c, ldc, UPLO, TR_c);

        if (UPLO == 'L') {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int j0 = jc + jb;
                if (TR_c == 'N') {
                    wgemm_(NN, CN, &trailing, &jb, &K, &alpha_c,
                           &A_(j0, 0), &lda, &A_(jc, 0), &lda,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(CN, NN, &trailing, &jb, &K, &alpha_c,
                           &A_(0, j0), &lda, &A_(0, jc), &lda,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                }
            }
        } else {
            if (jc > 0) {
                if (TR_c == 'N') {
                    wgemm_(NN, CN, &jc, &jb, &K, &alpha_c,
                           &A_(0, 0), &lda, &A_(jc, 0), &lda,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(CN, NN, &jc, &jb, &K, &alpha_c,
                           &A_(0, 0), &lda, &A_(0, jc), &lda,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef C_
