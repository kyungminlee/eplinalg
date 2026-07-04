/*
 * Fortran-callable helpers that drive the C bridge. They exist solely
 * for the parity tests: a Fortran driver invokes
 * ``c_real_mumps_solve`` / ``c_complex_mumps_solve`` (resolved via
 * BIND(C, name=...)) to run a solve through the C bridge and compare
 * the result bit-for-bit against the same problem solved through the
 * Fortran-direct path.
 *
 * Both paths converge to the same migrated <prefix>mumps_f77_ Fortran
 * entry; if the C struct extraction in mumps_c.c diverges from the
 * Fortran-direct path, the parity test surfaces it as a non-zero diff.
 *
 * Per-target choices (header path + entry-point names + struct-type
 * names) are supplied by the build system through TARGET_REAL_HEADER /
 * TARGET_REAL_MUMPS_C / TARGET_REAL_STRUC_C macros (and their COMPLEX
 * counterparts) — see tests/mumps/CMakeLists.txt.
 */

#include TARGET_REAL_HEADER
#include TARGET_COMPLEX_HEADER
#include <mpi.h>

int c_real_mumps_solve(int *n_in, MUMPS_INT8 *nnz_in,
                       MUMPS_INT *irn, MUMPS_INT *jcn,
                       DMUMPS_REAL *a_vals, DMUMPS_REAL *rhs)
{
    TARGET_REAL_STRUC_C id = {0};
    int code, myid;

    /* Centralized input: the assembled matrix and RHS are supplied on the
     * host only; the centralized solution comes back in rhs on the host
     * only. MUMPS calls stay collective on every rank. */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    TARGET_REAL_MUMPS_C(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    if (myid == 0) {
        id.n   = *n_in;
        id.nnz = *nnz_in;
        id.irn = irn;
        id.jcn = jcn;
        id.a   = a_vals;
        id.rhs = rhs;
    }
    id.job = 6;
    TARGET_REAL_MUMPS_C(&id);
    code = id.infog[0];

    id.job = -2;
    TARGET_REAL_MUMPS_C(&id);
    return code;
}

int c_complex_mumps_solve(int *n_in, MUMPS_INT8 *nnz_in,
                          MUMPS_INT *irn, MUMPS_INT *jcn,
                          mumps_double_complex *a_vals,
                          mumps_double_complex *rhs)
{
    TARGET_COMPLEX_STRUC_C id = {0};
    int code, myid;

    /* Centralized input: matrix + RHS on the host only; centralized
     * solution returned in rhs on the host only. Calls stay collective. */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    id.par = 1;
    id.sym = 0;
    id.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    id.job = -1;
    TARGET_COMPLEX_MUMPS_C(&id);
    if (id.infog[0] < 0) return id.infog[0];

    id.icntl[0] = -1; id.icntl[1] = -1; id.icntl[2] = -1; id.icntl[3] = 0;

    if (myid == 0) {
        id.n   = *n_in;
        id.nnz = *nnz_in;
        id.irn = irn;
        id.jcn = jcn;
        id.a   = a_vals;
        id.rhs = rhs;
    }
    id.job = 6;
    TARGET_COMPLEX_MUMPS_C(&id);
    code = id.infog[0];

    id.job = -2;
    TARGET_COMPLEX_MUMPS_C(&id);
    return code;
}

