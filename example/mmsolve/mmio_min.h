/*
 * mmio_min — a minimal Matrix Market reader/writer for the mmsolve example.
 *
 * Deliberately small and self-contained (no dependency on NIST mmio). It
 * covers exactly what a sparse linear solve needs:
 *
 *   - sparse matrices: banner "matrix coordinate real|complex general|symmetric"
 *   - dense vectors:   banner "matrix array  real|complex general" (n x 1)
 *
 * Everything is stored internally as double (real part in `val`, imaginary in
 * `ival` when complex); the typed solve in mmsolve.c converts to the selected
 * MUMPS working precision (float/double for s/c/d/z, long double for e/y,
 * __float128 for q/x, double-double {d,0} promotion for m/w). Indices are kept
 * 1-based, exactly as Matrix Market stores them and as MUMPS expects.
 */
#ifndef MMIO_MIN_H
#define MMIO_MIN_H

typedef struct {
    int    is_complex;    /* 0 = real, 1 = complex                     */
    int    is_symmetric;  /* 1 = only lower triangle stored (sym)      */
    int    is_coordinate; /* 1 = sparse coordinate, 0 = dense array    */
    int    m, n;          /* rows, cols                                */
    long   nnz;           /* stored entries (== m for a dense n x 1)   */
    int   *irn, *jcn;     /* 1-based row/col indices (coordinate only) */
    double *val;          /* real parts  (length nnz)                  */
    double *ival;         /* imag parts  (length nnz) or NULL if real  */
} MM;

/* Read a Matrix Market file into `out`. Returns 0 on success, non-zero on
 * error (a message is printed to stderr). On success the caller owns the
 * arrays and must call mm_free(out). */
int mm_read(const char *path, MM *out);

/* Write a dense n x 1 vector as a Matrix Market "array" file. `im` may be
 * NULL for a real vector. Returns 0 on success. */
int mm_write_vector(const char *path, int n, const double *re, const double *im);

void mm_free(MM *m);

#endif /* MMIO_MIN_H */
