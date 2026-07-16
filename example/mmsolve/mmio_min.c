#include "mmio_min.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *skip_ws(char *s) { while (*s && isspace((unsigned char)*s)) s++; return s; }

/* Lower-case a token in place. */
static void lower(char *s) { for (; *s; s++) *s = (char)tolower((unsigned char)*s); }

/* Read the next non-comment, non-blank line into buf. Returns 0 on success. */
static int next_data_line(FILE *f, char *buf, size_t n)
{
    while (fgets(buf, (int)n, f)) {
        char *p = skip_ws(buf);
        if (*p == '\0' || *p == '%') continue;  /* blank or comment */
        return 0;
    }
    return -1;
}

int mm_read(const char *path, MM *out)
{
    FILE *f = fopen(path, "r");
    char line[8192];
    if (!f) { fprintf(stderr, "mm_read: cannot open %s\n", path); return 1; }

    memset(out, 0, sizeof *out);

    /* ── banner ──────────────────────────────────────────────────────
     * %%MatrixMarket matrix <coordinate|array> <real|complex|integer|pattern>
     *                       <general|symmetric|hermitian|skew-symmetric>     */
    if (!fgets(line, sizeof line, f)) { fprintf(stderr, "mm_read: empty file\n"); fclose(f); return 1; }
    {
        char tag[64] = "", object[64] = "", format[64] = "", field[64] = "", symm[64] = "";
        int got = sscanf(line, "%63s %63s %63s %63s %63s", tag, object, format, field, symm);
        lower(tag); lower(object); lower(format); lower(field); lower(symm);
        if (got < 4 || strcmp(tag, "%%matrixmarket") != 0 || strcmp(object, "matrix") != 0) {
            fprintf(stderr, "mm_read: %s: not a Matrix Market matrix file\n", path);
            fclose(f); return 1;
        }
        if      (strcmp(format, "coordinate") == 0) out->is_coordinate = 1;
        else if (strcmp(format, "array") == 0)      out->is_coordinate = 0;
        else { fprintf(stderr, "mm_read: %s: unsupported format '%s'\n", path, format); fclose(f); return 1; }

        if      (strcmp(field, "real") == 0 || strcmp(field, "integer") == 0) out->is_complex = 0;
        else if (strcmp(field, "complex") == 0)                                out->is_complex = 1;
        else { fprintf(stderr, "mm_read: %s: unsupported field '%s' (pattern not supported)\n", path, field); fclose(f); return 1; }

        /* Symmetry. MUMPS' sym=2 is (complex-)*symmetric* — A = A^T with only
         * the lower triangle stored — which maps cleanly onto Matrix Market's
         * "symmetric". It has no Hermitian (A = A^H) mode, so we reject
         * "hermitian" and "skew-symmetric" rather than silently mis-solve a
         * half-stored matrix. (Callers wanting those must expand to a full
         * "general" matrix themselves.) */
        const char *sy = (got >= 5) ? symm : "general";
        if      (strcmp(sy, "general") == 0)   out->is_symmetric = 0;
        else if (strcmp(sy, "symmetric") == 0) out->is_symmetric = 1;
        else {
            fprintf(stderr, "mm_read: %s: symmetry '%s' unsupported "
                    "(use 'general', or 'symmetric' for MUMPS sym=2)\n", path, sy);
            fclose(f); return 1;
        }
    }

    /* ── size line ───────────────────────────────────────────────────*/
    if (next_data_line(f, line, sizeof line)) { fprintf(stderr, "mm_read: %s: missing size line\n", path); fclose(f); return 1; }

    if (out->is_coordinate) {
        long nnz = 0;
        if (sscanf(line, "%d %d %ld", &out->m, &out->n, &nnz) != 3) {
            fprintf(stderr, "mm_read: %s: bad coordinate size line\n", path); fclose(f); return 1;
        }
        out->nnz  = nnz;
        out->irn  = malloc(sizeof(int) * (size_t)nnz);
        out->jcn  = malloc(sizeof(int) * (size_t)nnz);
        out->val  = malloc(sizeof(double) * (size_t)nnz);
        out->ival = out->is_complex ? malloc(sizeof(double) * (size_t)nnz) : NULL;
        if (!out->irn || !out->jcn || !out->val || (out->is_complex && !out->ival)) {
            fprintf(stderr, "mm_read: out of memory\n"); mm_free(out); fclose(f); return 1;
        }
        for (long k = 0; k < nnz; k++) {
            if (next_data_line(f, line, sizeof line)) {
                fprintf(stderr, "mm_read: %s: expected %ld entries, got %ld\n", path, nnz, k);
                mm_free(out); fclose(f); return 1;
            }
            int i, j; double re = 0.0, im = 0.0;
            if (out->is_complex) {
                if (sscanf(line, "%d %d %lf %lf", &i, &j, &re, &im) != 4) {
                    fprintf(stderr, "mm_read: %s: bad complex entry at line %ld\n", path, k); mm_free(out); fclose(f); return 1;
                }
                out->ival[k] = im;
            } else {
                if (sscanf(line, "%d %d %lf", &i, &j, &re) != 3) {
                    fprintf(stderr, "mm_read: %s: bad real entry at line %ld\n", path, k); mm_free(out); fclose(f); return 1;
                }
            }
            out->irn[k] = i; out->jcn[k] = j; out->val[k] = re;
        }
    } else {
        /* dense array: "M N" then M*N values column-major. We only need N==1. */
        int cols = 1;
        if (sscanf(line, "%d %d", &out->m, &cols) < 1) {
            fprintf(stderr, "mm_read: %s: bad array size line\n", path); fclose(f); return 1;
        }
        out->n = cols;
        if (cols != 1) { fprintf(stderr, "mm_read: %s: only n x 1 vectors supported (got %d cols)\n", path, cols); fclose(f); return 1; }
        long cnt = out->m;
        out->nnz  = cnt;
        out->val  = malloc(sizeof(double) * (size_t)cnt);
        out->ival = out->is_complex ? malloc(sizeof(double) * (size_t)cnt) : NULL;
        if (!out->val || (out->is_complex && !out->ival)) {
            fprintf(stderr, "mm_read: out of memory\n"); mm_free(out); fclose(f); return 1;
        }
        for (long k = 0; k < cnt; k++) {
            if (next_data_line(f, line, sizeof line)) {
                fprintf(stderr, "mm_read: %s: expected %ld values, got %ld\n", path, cnt, k); mm_free(out); fclose(f); return 1;
            }
            double re = 0.0, im = 0.0;
            if (out->is_complex) {
                if (sscanf(line, "%lf %lf", &re, &im) != 2) { fprintf(stderr, "mm_read: %s: bad complex value at %ld\n", path, k); mm_free(out); fclose(f); return 1; }
                out->ival[k] = im;
            } else {
                if (sscanf(line, "%lf", &re) != 1) { fprintf(stderr, "mm_read: %s: bad value at %ld\n", path, k); mm_free(out); fclose(f); return 1; }
            }
            out->val[k] = re;
        }
    }

    fclose(f);
    return 0;
}

int mm_write_vector(const char *path, int n, const double *re, const double *im)
{
    FILE *f = fopen(path, "w");
    if (!f) { fprintf(stderr, "mm_write_vector: cannot open %s\n", path); return 1; }
    fprintf(f, "%%%%MatrixMarket matrix array %s general\n", im ? "complex" : "real");
    fprintf(f, "%% solution vector written by mmsolve\n");
    fprintf(f, "%d 1\n", n);
    for (int i = 0; i < n; i++) {
        if (im) fprintf(f, "%.17g %.17g\n", re[i], im[i]);
        else    fprintf(f, "%.17g\n", re[i]);
    }
    fclose(f);
    return 0;
}

void mm_free(MM *m)
{
    if (!m) return;
    free(m->irn); free(m->jcn); free(m->val); free(m->ival);
    m->irn = m->jcn = NULL; m->val = m->ival = NULL;
}
