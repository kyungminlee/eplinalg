/*
 * wsyr2k — multifloats complex (DD) symmetric rank-2k.
 *   C := alpha · (A · Bᵀ + B · Aᵀ) + beta · C        (TRANS='N')
 *   C := alpha · (Aᵀ · B + Bᵀ · A) + beta · C        (TRANS='T'/'C')
 * Blocked: scalar diagonal + two wgemm trailing calls per off-diag.
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

#define WSYR2K_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int syr2k_nb(void) { if (g_nb == 0) g_nb = env_int("WSYR2K_NB", 64); return g_nb; }

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
#define B_(i, j)  b[static_cast<std::size_t>(j) * ldb + (i)]
#define C_(i, j)  c[static_cast<std::size_t>(j) * ldc + (i)]

void syr2k_diag_add(int jc, int jb, int K, T alpha,
                    const T *a, int lda,
                    const T *b, int ldb,
                    T *c, int ldc,
                    char UPLO, char TR)
{
    if (TR == 'N') {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            for (int l = 0; l < K; ++l) {
                const T t1 = cmul(alpha, A_(j, l));
                const T t2 = cmul(alpha, B_(j, l));
                for (int i = i_lo; i < i_hi; ++i) {
                    cj[i] = cadd(cj[i], cadd(cmul(B_(i, l), t1), cmul(A_(i, l), t2)));
                }
            }
        }
    } else {
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            const T *Aj = a + static_cast<std::size_t>(j) * lda;
            const T *Bj = b + static_cast<std::size_t>(j) * ldb;
            for (int i = i_lo; i < i_hi; ++i) {
                const T *Ai = a + static_cast<std::size_t>(i) * lda;
                const T *Bi = b + static_cast<std::size_t>(i) * ldb;
                T s = zero_cdd;
                for (int l = 0; l < K; ++l) {
                    s = cadd(s, cadd(cmul(Ai[l], Bj[l]), cmul(Bi[l], Aj[l])));
                }
                cj[i] = cadd(cj[i], cmul(alpha, s));
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void wsyr2k_(
    const char *uplo, const char *trans,
    const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    std::size_t uplo_len, std::size_t trans_len)
{
    (void)uplo_len; (void)trans_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const char UPLO = up(uplo);
    char TR = up(trans);
    if (TR == 'C') TR = 'T';

    if (N == 0) return;

    const char NN[1] = {'N'};
    const char TN[1] = {'T'};

    if (cdd_iszero(alpha) || K == 0) {
        if (cdd_isone(beta)) return;
#ifdef _OPENMP
        const bool use_omp = (N >= WSYR2K_OMP_MIN && omp_get_max_threads() > 1);
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

    const int nb = syr2k_nb();

#ifdef _OPENMP
    const bool use_omp = (N >= WSYR2K_OMP_MIN && omp_get_max_threads() > 1);
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

        syr2k_diag_add(jc, jb, K, alpha, a, lda, b, ldb, c, ldc, UPLO, TR);

        if (UPLO == 'L') {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int j0 = jc + jb;
                if (TR == 'N') {
                    wgemm_(NN, TN, &trailing, &jb, &K, &alpha,
                           &A_(j0, 0), &lda, &B_(jc, 0), &ldb,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                    wgemm_(NN, TN, &trailing, &jb, &K, &alpha,
                           &B_(j0, 0), &ldb, &A_(jc, 0), &lda,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(TN, NN, &trailing, &jb, &K, &alpha,
                           &A_(0, j0), &lda, &B_(0, jc), &ldb,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                    wgemm_(TN, NN, &trailing, &jb, &K, &alpha,
                           &B_(0, j0), &ldb, &A_(0, jc), &lda,
                           &one_cdd, &C_(j0, jc), &ldc, 1, 1);
                }
            }
        } else {
            if (jc > 0) {
                if (TR == 'N') {
                    wgemm_(NN, TN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &B_(jc, 0), &ldb,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                    wgemm_(NN, TN, &jc, &jb, &K, &alpha,
                           &B_(0, 0), &ldb, &A_(jc, 0), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(TN, NN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &B_(0, jc), &ldb,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                    wgemm_(TN, NN, &jc, &jb, &K, &alpha,
                           &B_(0, 0), &ldb, &A_(0, jc), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
