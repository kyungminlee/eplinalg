/*
 * wsyrk — multifloats complex (DD) symmetric rank-k update. TRANS ∈ {N, T}.
 * Blocked: scalar diagonal + wgemm trailing.
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

#define WSYRK_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int syrk_nb(void) { if (g_nb == 0) g_nb = env_int("WSYRK_NB", 64); return g_nb; }

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

void syrk_diag_add(int jc, int jb, int K, T alpha,
                   const T *a, int lda, T *c, int ldc,
                   char UPLO, char TR)
{
    if (TR == 'N') {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            for (int l = 0; l < K; ++l) {
                const T ajl = A_(j, l);
                if (!cdd_iszero(ajl)) {
                    const T t = cmul(alpha, ajl);
                    for (int i = i_lo; i < i_hi; ++i) cj[i] = cadd(cj[i], cmul(t, A_(i, l)));
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
                T s = zero_cdd;
                for (int l = 0; l < K; ++l) s = cadd(s, cmul(Ai[l], Aj[l]));
                cj[i] = cadd(cj[i], cmul(alpha, s));
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void wsyrk_(
    const char *uplo, const char *trans,
    const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *beta_,
    T *c, const int *ldc_,
    std::size_t uplo_len, std::size_t trans_len)
{
    (void)uplo_len; (void)trans_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);
    const char TR = up(trans);

    if (N == 0) return;

    const char NN[1] = {'N'};
    const char TN[1] = {'T'};

    if (cdd_iszero(alpha) || K == 0) {
        if (cdd_isone(beta)) return;
#ifdef _OPENMP
        const bool use_omp = (N >= WSYRK_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (cdd_iszero(beta)) for (int i = i_lo; i < i_hi; ++i) cj[i] = zero_cdd;
            else                  for (int i = i_lo; i < i_hi; ++i) cj[i] = cmul(cj[i], beta);
        }
        return;
    }

    const int nb = syrk_nb();

#ifdef _OPENMP
    const bool use_omp = (N >= WSYRK_OMP_MIN && omp_get_max_threads() > 1);
    #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
    for (int jc = 0; jc < N; jc += nb) {
        const int jb = (N - jc < nb) ? (N - jc) : nb;

        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j : 0;
            const int i_hi = (UPLO == 'L') ? N : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (cdd_iszero(beta))      for (int i = i_lo; i < i_hi; ++i) cj[i] = zero_cdd;
            else if (!cdd_isone(beta)) for (int i = i_lo; i < i_hi; ++i) cj[i] = cmul(cj[i], beta);
        }

        syrk_diag_add(jc, jb, K, alpha, a, lda, c, ldc, UPLO, TR);

        if (UPLO == 'L') {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int j0 = jc + jb;
                if (TR == 'N') {
                    wgemm_(NN, TN, &trailing, &jb, &K, &alpha,
                           &A_(j0, 0), &lda, &A_(jc, 0), &lda,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(TN, NN, &trailing, &jb, &K, &alpha,
                           &A_(0, j0), &lda, &A_(0, jc), &lda,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                }
            }
        } else {
            if (jc > 0) {
                if (TR == 'N') {
                    wgemm_(NN, TN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &A_(jc, 0), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(TN, NN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &A_(0, jc), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef C_
