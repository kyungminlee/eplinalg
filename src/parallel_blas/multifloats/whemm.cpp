/*
 * whemm — multifloats complex (DD) Hermitian matrix multiply.
 * Blocked + wgemm trailing. Reflection gemms use 'C' (conjugate
 * transpose); scalar diagonal kernel keeps the diagonal real.
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

#define WHEMM_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int hemm_nb(void) { if (g_nb == 0) g_nb = env_int("WHEMM_NB", 64); return g_nb; }

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const R rzero{0.0, 0.0};
const T zero_cdd{ rzero, rzero };
const T one_cdd { R{1.0, 0.0}, rzero };

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

void hemm_diag_add_L(int ic, int ib, int N, T alpha,
                     const T *a, int lda, const T *b, int ldb,
                     T *c, int ldc, char UPLO)
{
    for (int j = 0; j < N; ++j) {
        T *cj = c + static_cast<std::size_t>(j) * ldc;
        const T *bj = b + static_cast<std::size_t>(j) * ldb;
        if (UPLO == 'L') {
            for (int i = ic; i < ic + ib; ++i) {
                const T temp1 = cmul(alpha, bj[i]);
                T temp2 = zero_cdd;
                for (int k = ic; k < i; ++k) {
                    cj[k]  = cadd(cj[k], cmul(temp1, cconj(A_(i, k))));
                    temp2  = cadd(temp2, cmul(bj[k], A_(i, k)));
                }
                /* Diagonal of A is real: take real part. */
                const T diag = T{ A_(i, i).re, rzero };
                cj[i] = cadd(cj[i], cadd(cmul(temp1, diag), cmul(alpha, temp2)));
            }
        } else {
            for (int i = ic + ib - 1; i >= ic; --i) {
                const T temp1 = cmul(alpha, bj[i]);
                T temp2 = zero_cdd;
                for (int k = i + 1; k < ic + ib; ++k) {
                    cj[k]  = cadd(cj[k], cmul(temp1, cconj(A_(i, k))));
                    temp2  = cadd(temp2, cmul(bj[k], A_(i, k)));
                }
                const T diag = T{ A_(i, i).re, rzero };
                cj[i] = cadd(cj[i], cadd(cmul(temp1, diag), cmul(alpha, temp2)));
            }
        }
    }
}

void hemm_diag_add_R(int jc, int jb, int M, T alpha,
                     const T *a, int lda, const T *b, int ldb,
                     T *c, int ldc, char UPLO)
{
    for (int j = jc; j < jc + jb; ++j) {
        T *cj = c + static_cast<std::size_t>(j) * ldc;
        {
            const T diag = T{ A_(j, j).re, rzero };
            const T t = cmul(alpha, diag);
            for (int i = 0; i < M; ++i) cj[i] = cadd(cj[i], cmul(t, B_(i, j)));
        }
        if (UPLO == 'L') {
            for (int k = jc; k < j; ++k) {
                const T t = cmul(alpha, cconj(A_(j, k)));
                if (!cdd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cadd(cj[i], cmul(t, B_(i, k)));
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = cmul(alpha, A_(k, j));
                if (!cdd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cadd(cj[i], cmul(t, B_(i, k)));
            }
        } else {
            for (int k = jc; k < j; ++k) {
                const T t = cmul(alpha, A_(k, j));
                if (!cdd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cadd(cj[i], cmul(t, B_(i, k)));
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = cmul(alpha, cconj(A_(j, k)));
                if (!cdd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cadd(cj[i], cmul(t, B_(i, k)));
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void whemm_(
    const char *side, const char *uplo,
    const int *m_, const int *n_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    std::size_t side_len, std::size_t uplo_len)
{
    (void)side_len; (void)uplo_len;
    const int M = *m_, N = *n_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const char SIDE = up(side);
    const char UPLO = up(uplo);

    if (M == 0 || N == 0) return;

    const char NN[1] = {'N'};
    const char CN[1] = {'C'};

    if (cdd_iszero(alpha)) {
        if (cdd_isone(beta)) return;
#ifdef _OPENMP
        const int axis = (SIDE == 'L') ? M : N;
        const bool use_omp = (axis >= WHEMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (cdd_iszero(beta)) for (int i = 0; i < M; ++i) cj[i] = zero_cdd;
            else                  for (int i = 0; i < M; ++i) cj[i] = cmul(cj[i], beta);
        }
        return;
    }

    const int nb = hemm_nb();

    if (SIDE == 'L') {
#ifdef _OPENMP
        const bool use_omp = (M >= WHEMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            for (int j = 0; j < N; ++j) {
                T *cj = c + static_cast<std::size_t>(j) * ldc;
                if (cdd_iszero(beta))      for (int i = ic; i < ic + ib; ++i) cj[i] = zero_cdd;
                else if (!cdd_isone(beta)) for (int i = ic; i < ic + ib; ++i) cj[i] = cmul(cj[i], beta);
            }
            if (UPLO == 'L') {
                if (ic > 0) {
                    wgemm_(NN, NN, &ib, &N, &ic, &alpha,
                           &A_(ic, 0), &lda, &B_(0, 0), &ldb,
                           &one_cdd, &C_(ic, 0), &ldc, 1, 1);
                }
                hemm_diag_add_L(ic, ib, N, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = M - ic - ib;
                if (trailing > 0) {
                    wgemm_(CN, NN, &ib, &N, &trailing, &alpha,
                           &A_(ic + ib, ic), &lda, &B_(ic + ib, 0), &ldb,
                           &one_cdd, &C_(ic, 0), &ldc, 1, 1);
                }
            } else {
                if (ic > 0) {
                    wgemm_(CN, NN, &ib, &N, &ic, &alpha,
                           &A_(0, ic), &lda, &B_(0, 0), &ldb,
                           &one_cdd, &C_(ic, 0), &ldc, 1, 1);
                }
                hemm_diag_add_L(ic, ib, N, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = M - ic - ib;
                if (trailing > 0) {
                    wgemm_(NN, NN, &ib, &N, &trailing, &alpha,
                           &A_(ic, ic + ib), &lda, &B_(ic + ib, 0), &ldb,
                           &one_cdd, &C_(ic, 0), &ldc, 1, 1);
                }
            }
        }
    } else {
#ifdef _OPENMP
        const bool use_omp = (N >= WHEMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            for (int j = jc; j < jc + jb; ++j) {
                T *cj = c + static_cast<std::size_t>(j) * ldc;
                if (cdd_iszero(beta))      for (int i = 0; i < M; ++i) cj[i] = zero_cdd;
                else if (!cdd_isone(beta)) for (int i = 0; i < M; ++i) cj[i] = cmul(cj[i], beta);
            }
            if (UPLO == 'L') {
                if (jc > 0) {
                    wgemm_(NN, CN, &M, &jb, &jc, &alpha,
                           &B_(0, 0), &ldb, &A_(jc, 0), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
                hemm_diag_add_R(jc, jb, M, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = N - jc - jb;
                if (trailing > 0) {
                    wgemm_(NN, NN, &M, &jb, &trailing, &alpha,
                           &B_(0, jc + jb), &ldb, &A_(jc + jb, jc), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
            } else {
                if (jc > 0) {
                    wgemm_(NN, NN, &M, &jb, &jc, &alpha,
                           &B_(0, 0), &ldb, &A_(0, jc), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
                hemm_diag_add_R(jc, jb, M, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = N - jc - jb;
                if (trailing > 0) {
                    wgemm_(NN, CN, &M, &jb, &trailing, &alpha,
                           &B_(0, jc + jb), &ldb, &A_(jc, jc + jb), &lda,
                           &one_cdd, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
