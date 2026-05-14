/*
 * egemm — kind10 (REAL(KIND=10), x86-64 80-bit long double) GEMM overlay.
 *
 * Step 0: naive triple-loop reference implementation. Purpose is to
 * (a) prove the Fortran-ABI C boundary end-to-end and (b) let the
 * existing tests/blas/level3/test_egemm exercise the overlay path
 * from day one. This is the correctness floor — every subsequent
 * step (blocking, packing, OpenMP) must remain consistent with this
 * implementation modulo summation order, within the 10-ulp tolerance.
 *
 * Real Goto-style blocked + packed implementation lands in step 3
 * (serial) and step 4 (OpenMP-parallel), replacing this body.
 *
 * Fortran ABI (gfortran):
 *   - subroutine name lowercased + trailing underscore: `egemm_`
 *   - scalars passed by pointer
 *   - character args followed by hidden trailing `size_t` lengths
 *   - REAL(KIND=10) ↔ long double (x86-64 80-bit extended)
 *
 * Column-major storage:
 *   A(i,j) at index i + j*lda (0-based) / 1 + i + (j-1)*lda (1-based)
 */

#include <stddef.h>
#include <ctype.h>

typedef long double T;

static int trans_code(const char *p, size_t len)
{
    /* gfortran passes hidden length; we only look at the first char. */
    (void)len;
    char c = (char)toupper((unsigned char)*p);
    return c;  /* 'N' | 'T' | 'C' (C == T for real types) */
}

void egemm_(
    const char *transa, const char *transb,
    const int *m_, const int *n_, const int *k_,
    const T *alpha_,
    const T *a, const int *lda_,
    const T *b, const int *ldb_,
    const T *beta_,
    T *c, const int *ldc_,
    size_t transa_len, size_t transb_len)
{
    const int m = *m_, n = *n_, k = *k_;
    const int lda = *lda_, ldb = *ldb_, ldc = *ldc_;
    const T alpha = *alpha_, beta = *beta_;
    const int ta = trans_code(transa, transa_len);
    const int tb = trans_code(transb, transb_len);

    /* C <- beta * C (the alpha=0 quick-return path). */
    if (alpha == 0.0L || k == 0) {
        for (int j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            if (beta == 0.0L) {
                for (int i = 0; i < m; ++i) cj[i] = 0.0L;
            } else if (beta != 1.0L) {
                for (int i = 0; i < m; ++i) cj[i] *= beta;
            }
        }
        return;
    }

    /* Apply beta to C up front; then accumulate alpha*op(A)*op(B). */
    for (int j = 0; j < n; ++j) {
        T *cj = &c[(size_t)j * ldc];
        if (beta == 0.0L) {
            for (int i = 0; i < m; ++i) cj[i] = 0.0L;
        } else if (beta != 1.0L) {
            for (int i = 0; i < m; ++i) cj[i] *= beta;
        }
    }

    /* Naive ijk loops, specialized on transpose flags. */
    if (ta == 'N' && tb == 'N') {
        for (int j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            const T *bj = &b[(size_t)j * ldb];
            for (int p = 0; p < k; ++p) {
                const T t = alpha * bj[p];
                const T *ap = &a[(size_t)p * lda];
                for (int i = 0; i < m; ++i) cj[i] += t * ap[i];
            }
        }
    } else if (ta == 'N' && tb != 'N') {
        for (int j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            for (int p = 0; p < k; ++p) {
                const T t = alpha * b[(size_t)p * ldb + j];
                const T *ap = &a[(size_t)p * lda];
                for (int i = 0; i < m; ++i) cj[i] += t * ap[i];
            }
        }
    } else if (ta != 'N' && tb == 'N') {
        for (int j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            const T *bj = &b[(size_t)j * ldb];
            for (int i = 0; i < m; ++i) {
                const T *ai = &a[(size_t)i * lda];
                T acc = 0.0L;
                for (int p = 0; p < k; ++p) acc += ai[p] * bj[p];
                cj[i] += alpha * acc;
            }
        }
    } else {
        for (int j = 0; j < n; ++j) {
            T *cj = &c[(size_t)j * ldc];
            for (int i = 0; i < m; ++i) {
                const T *ai = &a[(size_t)i * lda];
                T acc = 0.0L;
                for (int p = 0; p < k; ++p)
                    acc += ai[p] * b[(size_t)p * ldb + j];
                cj[i] += alpha * acc;
            }
        }
    }
}
