/*
 * wher2k — multifloats complex (DD) Hermitian rank-2k.
 *   C := alpha · A · Bᴴ + conj(alpha) · B · Aᴴ + beta · C  (TRANS='N')
 *   C := alpha · Aᴴ · B + conj(alpha) · Bᴴ · A + beta · C  (TRANS='C')
 * alpha complex, beta real. Diagonal of C stays real.
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

#define WHER2K_OMP_MIN 32

int env_int(const char *name, int dflt) {
    const char *s = std::getenv(name);
    if (!s || !*s) return dflt;
    int v = std::atoi(s);
    return v > 0 ? v : dflt;
}
int g_nb = 0;
int her2k_nb(void) { if (g_nb == 0) g_nb = env_int("WHER2K_NB", 64); return g_nb; }

inline char up(const char *p) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(*p)));
}

const R rzero{0.0, 0.0};
const R rone {1.0, 0.0};
const T czero{ rzero, rzero };
const T cone { rone,  rzero };

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
#define B_(i, j)  b[static_cast<std::size_t>(j) * ldb + (i)]
#define C_(i, j)  c[static_cast<std::size_t>(j) * ldc + (i)]

void her2k_diag_add(int jc, int jb, int K, T alpha,
                    const T *a, int lda,
                    const T *b, int ldb,
                    T *c, int ldc,
                    char UPLO, char TR_c)
{
    const T alpha_conj = cconj(alpha);
    if (TR_c == 'N') {
        /* C(I,J) += α A(I,l) conj(B(J,l)) + conj(α) B(I,l) conj(A(J,l)) */
        for (int j = jc; j < jc + jb; ++j) {
            const int i_lo = (UPLO == 'L') ? j     : jc;
            const int i_hi = (UPLO == 'L') ? jc+jb : j + 1;
            T *cj = c + static_cast<std::size_t>(j) * ldc;
            for (int l = 0; l < K; ++l) {
                const T t1 = cmul(alpha,       cconj(B_(j, l)));
                const T t2 = cmul(alpha_conj,  cconj(A_(j, l)));
                for (int i = i_lo; i < i_hi; ++i) {
                    const T prod = cadd(cmul(A_(i, l), t1), cmul(B_(i, l), t2));
                    if (i == j) cj[i] = T{ cj[i].re + prod.re, cj[i].im };
                    else        cj[i] = cadd(cj[i], prod);
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
                T s1 = czero, s2 = czero;
                for (int l = 0; l < K; ++l) {
                    s1 = cadd(s1, cmul(cconj(Ai[l]), Bj[l]));
                    s2 = cadd(s2, cmul(cconj(Bi[l]), Aj[l]));
                }
                const T as = cadd(cmul(alpha, s1), cmul(alpha_conj, s2));
                if (i == j) cj[i] = T{ cj[i].re + as.re, cj[i].im };
                else        cj[i] = cadd(cj[i], as);
            }
        }
    }
}

} /* anonymous namespace */

extern "C" void wher2k_(
    const char *uplo, const char *trans,
    const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const R *beta_,
    T *c, const int *ldc_,
    std::size_t uplo_len, std::size_t trans_len)
{
    (void)uplo_len; (void)trans_len;
    const int N = *n_, K = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_;
    const R beta  = *beta_;
    const char UPLO = up(uplo);
    const char TR_c = up(trans);

    if (N == 0) return;

    const char NN[1] = {'N'};
    const char CN[1] = {'C'};
    const T alpha_conj = cconj(alpha);

    if (cdd_iszero(alpha) || K == 0) {
        if (dd_isone(beta)) {
            for (int j = 0; j < N; ++j) c[static_cast<std::size_t>(j) * ldc + j].im = rzero;
            return;
        }
#ifdef _OPENMP
        const bool use_omp = (N >= WHER2K_OMP_MIN && omp_get_max_threads() > 1);
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

    const int nb = her2k_nb();

#ifdef _OPENMP
    const bool use_omp = (N >= WHER2K_OMP_MIN && omp_get_max_threads() > 1);
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

        her2k_diag_add(jc, jb, K, alpha, a, lda, b, ldb, c, ldc, UPLO, TR_c);

        if (UPLO == 'L') {
            const int trailing = N - jc - jb;
            if (trailing > 0) {
                const int j0 = jc + jb;
                if (TR_c == 'N') {
                    wgemm_(NN, CN, &trailing, &jb, &K, &alpha,
                           &A_(j0, 0), &lda, &B_(jc, 0), &ldb,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                    wgemm_(NN, CN, &trailing, &jb, &K, &alpha_conj,
                           &B_(j0, 0), &ldb, &A_(jc, 0), &lda,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(CN, NN, &trailing, &jb, &K, &alpha,
                           &A_(0, j0), &lda, &B_(0, jc), &ldb,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                    wgemm_(CN, NN, &trailing, &jb, &K, &alpha_conj,
                           &B_(0, j0), &ldb, &A_(0, jc), &lda,
                           &cone, &C_(j0, jc), &ldc, 1, 1);
                }
            }
        } else {
            if (jc > 0) {
                if (TR_c == 'N') {
                    wgemm_(NN, CN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &B_(jc, 0), &ldb,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                    wgemm_(NN, CN, &jc, &jb, &K, &alpha_conj,
                           &B_(0, 0), &ldb, &A_(jc, 0), &lda,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                } else {
                    wgemm_(CN, NN, &jc, &jb, &K, &alpha,
                           &A_(0, 0), &lda, &B_(0, jc), &ldb,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                    wgemm_(CN, NN, &jc, &jb, &K, &alpha_conj,
                           &B_(0, 0), &ldb, &A_(0, jc), &lda,
                           &cone, &C_(0, jc), &ldc, 1, 1);
                }
            }
        }
    }
}

#undef A_
#undef B_
#undef C_
