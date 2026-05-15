/*
 * msymm — multifloats real (DD) symmetric matrix multiply.
 * Blocked: scalar diagonal + mgemm trailing.
 */

#include <cstddef>
#include <cstdlib>
#include <cctype>
#include <multifloats.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace mf = multifloats;
using T = mf::float64x2;

namespace {

#define MSYMM_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int symm_nb(void) { if (g_nb == 0) g_nb = env_int("MSYMM_NB", 64); return g_nb; }

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const T zero_dd{0.0, 0.0};
const T one_dd {1.0, 0.0};
inline bool dd_iszero(T x) { return x.limbs[0] == 0.0 && x.limbs[1] == 0.0; }
inline bool dd_isone (T x) { return x.limbs[0] == 1.0 && x.limbs[1] == 0.0; }

extern "C" void mgemm_(
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

void symm_diag_add_L(int ic, int ib, int N, T alpha,
                     const T *a, int lda, const T *b, int ldb,
                     T *c, int ldc, char UPLO)
{
    for (int j = 0; j < N; ++j) {
        T *cj = c + static_cast<std::size_t>(j) * ldc;
        const T *bj = b + static_cast<std::size_t>(j) * ldb;
        if (UPLO == 'L') {
            for (int i = ic; i < ic + ib; ++i) {
                const T temp1 = alpha * bj[i];
                T temp2 = zero_dd;
                for (int k = ic; k < i; ++k) {
                    cj[k]  = cj[k]  + temp1 * A_(i, k);
                    temp2  = temp2  + bj[k] * A_(i, k);
                }
                cj[i] = cj[i] + temp1 * A_(i, i) + alpha * temp2;
            }
        } else {
            for (int i = ic + ib - 1; i >= ic; --i) {
                const T temp1 = alpha * bj[i];
                T temp2 = zero_dd;
                for (int k = i + 1; k < ic + ib; ++k) {
                    cj[k]  = cj[k]  + temp1 * A_(i, k);
                    temp2  = temp2  + bj[k] * A_(i, k);
                }
                cj[i] = cj[i] + temp1 * A_(i, i) + alpha * temp2;
            }
        }
    }
}

void symm_diag_add_R(int jc, int jb, int M, T alpha,
                     const T *a, int lda, const T *b, int ldb,
                     T *c, int ldc, char UPLO)
{
    for (int j = jc; j < jc + jb; ++j) {
        T *cj = c + static_cast<std::size_t>(j) * ldc;
        {
            const T t = alpha * A_(j, j);
            for (int i = 0; i < M; ++i) cj[i] = cj[i] + t * B_(i, j);
        }
        if (UPLO == 'L') {
            for (int k = jc; k < j; ++k) {
                const T t = alpha * A_(j, k);
                if (!dd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cj[i] + t * B_(i, k);
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = alpha * A_(k, j);
                if (!dd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cj[i] + t * B_(i, k);
            }
        } else {
            for (int k = jc; k < j; ++k) {
                const T t = alpha * A_(k, j);
                if (!dd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cj[i] + t * B_(i, k);
            }
            for (int k = j + 1; k < jc + jb; ++k) {
                const T t = alpha * A_(j, k);
                if (!dd_iszero(t)) for (int i = 0; i < M; ++i) cj[i] = cj[i] + t * B_(i, k);
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void msymm_(
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
    const char TN[1] = {'T'};

    if (dd_iszero(alpha)) {
        if (dd_isone(beta)) return;
#ifdef _OPENMP
        const int axis = (SIDE == 'L') ? M : N;
        const bool use_omp = (axis >= MSYMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(static)
#endif
        for (int j = 0; j < N; ++j) {
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            if (dd_iszero(beta)) for (int i = 0; i < M; ++i) cj[i] = zero_dd;
            else                 for (int i = 0; i < M; ++i) cj[i] = cj[i] * beta;
        }
        return;
    }

    const int nb = symm_nb();

    if (SIDE == 'L') {
#ifdef _OPENMP
        const bool use_omp = (M >= MSYMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
        for (int ic = 0; ic < M; ic += nb) {
            const int ib = (M - ic < nb) ? (M - ic) : nb;
            for (int j = 0; j < N; ++j) {
                T *cj = c + static_cast<std::size_t>(j) * ldc;
                if (dd_iszero(beta))      for (int i = ic; i < ic + ib; ++i) cj[i] = zero_dd;
                else if (!dd_isone(beta)) for (int i = ic; i < ic + ib; ++i) cj[i] = cj[i] * beta;
            }
            if (UPLO == 'L') {
                if (ic > 0) {
                    mgemm_(NN, NN, &ib, &N, &ic, &alpha,
                           &A_(ic, 0), &lda, &B_(0, 0), &ldb,
                           &one_dd, &C_(ic, 0), &ldc, 1, 1);
                }
                symm_diag_add_L(ic, ib, N, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = M - ic - ib;
                if (trailing > 0) {
                    mgemm_(TN, NN, &ib, &N, &trailing, &alpha,
                           &A_(ic + ib, ic), &lda, &B_(ic + ib, 0), &ldb,
                           &one_dd, &C_(ic, 0), &ldc, 1, 1);
                }
            } else {
                if (ic > 0) {
                    mgemm_(TN, NN, &ib, &N, &ic, &alpha,
                           &A_(0, ic), &lda, &B_(0, 0), &ldb,
                           &one_dd, &C_(ic, 0), &ldc, 1, 1);
                }
                symm_diag_add_L(ic, ib, N, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = M - ic - ib;
                if (trailing > 0) {
                    mgemm_(NN, NN, &ib, &N, &trailing, &alpha,
                           &A_(ic, ic + ib), &lda, &B_(ic + ib, 0), &ldb,
                           &one_dd, &C_(ic, 0), &ldc, 1, 1);
                }
            }
        }
    } else {
#ifdef _OPENMP
        const bool use_omp = (N >= MSYMM_OMP_MIN && omp_get_max_threads() > 1);
        #pragma omp parallel for if(use_omp) schedule(dynamic, 1)
#endif
        for (int jc = 0; jc < N; jc += nb) {
            const int jb = (N - jc < nb) ? (N - jc) : nb;
            for (int j = jc; j < jc + jb; ++j) {
                T *cj = c + static_cast<std::size_t>(j) * ldc;
                if (dd_iszero(beta))      for (int i = 0; i < M; ++i) cj[i] = zero_dd;
                else if (!dd_isone(beta)) for (int i = 0; i < M; ++i) cj[i] = cj[i] * beta;
            }
            if (UPLO == 'L') {
                if (jc > 0) {
                    mgemm_(NN, TN, &M, &jb, &jc, &alpha,
                           &B_(0, 0), &ldb, &A_(jc, 0), &lda,
                           &one_dd, &C_(0, jc), &ldc, 1, 1);
                }
                symm_diag_add_R(jc, jb, M, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = N - jc - jb;
                if (trailing > 0) {
                    mgemm_(NN, NN, &M, &jb, &trailing, &alpha,
                           &B_(0, jc + jb), &ldb, &A_(jc + jb, jc), &lda,
                           &one_dd, &C_(0, jc), &ldc, 1, 1);
                }
            } else {
                if (jc > 0) {
                    mgemm_(NN, NN, &M, &jb, &jc, &alpha,
                           &B_(0, 0), &ldb, &A_(0, jc), &lda,
                           &one_dd, &C_(0, jc), &ldc, 1, 1);
                }
                symm_diag_add_R(jc, jb, M, alpha, a, lda, b, ldb, c, ldc, UPLO);
                const int trailing = N - jc - jb;
                if (trailing > 0) {
                    mgemm_(NN, TN, &M, &jb, &trailing, &alpha,
                           &B_(0, jc + jb), &ldb, &A_(jc, jc + jb), &lda,
                           &one_dd, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
