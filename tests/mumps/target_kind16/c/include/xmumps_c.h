/*
 * Self-contained MUMPS C header for the migrated xmumps archive
 * (``__float128``; complex arithmetic; prefix
 * x/X, maps to upstream ``z``).
 *
 * Declares XMUMPS_STRUC_C and xmumps_c() directly instead of
 * macro-renaming and #include-ing upstream ``zmumps_c.h``. The struct
 * layout is copied verbatim from MUMPS 5.8.2 ``zmumps_c.h``; its field
 * widths use this target's XMUMPS_REAL / XMUMPS_COMPLEX (mumps_float128 (scalar) / mumps_float128_complex (complex)),
 * defined below. This makes the struct ABI-identical to the xmumps
 * bridge object (upstream ``mumps_c.c`` compiled with
 * -DZMUMPS_STRUC_C=XMUMPS_STRUC_C -Dzmumps_c=xmumps_c and force-including
 * ``mumps_c_types_extended.h``, which overrides ZMUMPS_* to the same
 * width).
 *
 * Deliberately minimal: it defines only what the struct + prototype need
 * — MUMPS_INT / MUMPS_INT8 (mirroring upstream ``mumps_c_types.h`` via the
 * same ``mumps_int_def.h``), the widened widths, and MUMPS_CALL. It does
 * NOT pull in upstream ``mumps_c_types.h`` / ``mumps_c_types_extended.h``,
 * so a caller is not handed the MUMPS_ARITH_* flags, the S/CMUMPS_*
 * widths, or the mumps_complex / mumps_double_complex typedefs it never
 * uses.
 *
 * Derived from MUMPS 5.8.2 (Copyright 1991-2026 CERFACS, CNRS, ENS Lyon,
 * INP Toulouse, Inria, Mumps Technologies, University of Bordeaux),
 * released under the CeCILL-C license.
 */

#ifndef XMUMPS_C_H
#define XMUMPS_C_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mumps_compat.h"       /* MUMPS_CALL */

/* MUMPS_INT / MUMPS_INT8 — identical to upstream ``mumps_c_types.h`` (same
 * ``mumps_int_def.h``), without its S/C arithmetic widths, the
 * mumps_complex typedefs, or the MUMPS_ARITH_* selector flags. */
#include <stdint.h>
#include "mumps_int_def.h"      /* MUMPS_INTSIZE32 / MUMPS_INTSIZE64 */
#ifndef MUMPS_INT
# ifdef MUMPS_INTSIZE64
#  define MUMPS_INT int64_t
# else
#  define MUMPS_INT int
# endif
#endif
#ifndef MUMPS_INT8
# define MUMPS_INT8 int64_t
#endif

/* Extended real / complex widths (``__float128`` (REAL/COMPLEX(KIND=16))). */
#ifndef MUMPS_FLOAT128_TYPES
#define MUMPS_FLOAT128_TYPES
typedef __float128                  mumps_float128;
typedef struct { __float128 r, i; } mumps_float128_complex;
#endif
#define XMUMPS_REAL    mumps_float128
#define XMUMPS_COMPLEX mumps_float128_complex

#ifndef MUMPS_VERSION
/* Protected in case headers of other arithmetics are included */
#define MUMPS_VERSION "5.8.2"
#endif
#ifndef MUMPS_VERSION_MAX_LEN
#define MUMPS_VERSION_MAX_LEN 30
#endif

/*
 * Definition of the (simplified) MUMPS C structure.
 * NB: XMUMPS_COMPLEX are REAL types in s and d arithmetics.
 */
typedef struct {

    MUMPS_INT      sym, par, job;
    MUMPS_INT      comm_fortran;    /* Fortran communicator */
    MUMPS_INT      icntl[60];
    MUMPS_INT      keep[500];
    XMUMPS_REAL    cntl[15];
    XMUMPS_REAL    dkeep[230];
    MUMPS_INT8     keep8[150];
    MUMPS_INT      n;
    MUMPS_INT      nblk;

    MUMPS_INT      nz_alloc; /* used in matlab interface to decide if we
                                free + malloc when we have large variation */

    /* Assembled entry */
    MUMPS_INT      nz;
    MUMPS_INT8     nnz;
    MUMPS_INT      *irn;
    MUMPS_INT      *jcn;
    XMUMPS_COMPLEX *a;

    /* Distributed entry */
    MUMPS_INT      nz_loc;
    MUMPS_INT8     nnz_loc;
    MUMPS_INT      *irn_loc;
    MUMPS_INT      *jcn_loc;
    XMUMPS_COMPLEX *a_loc;

    /* Element entry */
    MUMPS_INT      nelt;
    MUMPS_INT      *eltptr;
    MUMPS_INT      *eltvar;
    XMUMPS_COMPLEX *a_elt;

    /* Matrix by blocks */
    MUMPS_INT      *blkptr;
    MUMPS_INT      *blkvar;

    /* Ordering, if given by user */
    MUMPS_INT      *perm_in;

    /* Orderings returned to user */
    MUMPS_INT      *sym_perm;    /* symmetric permutation */
    MUMPS_INT      *uns_perm;    /* column permutation */

    /* Scaling (inout but complicated) */
    XMUMPS_REAL    *colsca;
    XMUMPS_REAL    *rowsca;
    MUMPS_INT colsca_from_mumps;
    MUMPS_INT rowsca_from_mumps;

    /* Distributed scaling(out) */
    XMUMPS_REAL    *colsca_loc;
    XMUMPS_REAL    *rowsca_loc;

    /* Info after facto */
    MUMPS_INT      *rowind;
    MUMPS_INT      *colind;
    XMUMPS_COMPLEX *pivots;

    /* RHS, solution, ouptput data and statistics */
    XMUMPS_COMPLEX *rhs, *redrhs, *rhs_sparse, *sol_loc, *rhs_loc, *rhsintr;
    MUMPS_INT      *irhs_sparse, *irhs_ptr, *isol_loc, *irhs_loc, *glob2loc_rhs, *glob2loc_sol;
    MUMPS_INT      nrhs, lrhs, lredrhs, nz_rhs, lsol_loc, nloc_rhs, lrhs_loc, nsol_loc;
    MUMPS_INT      schur_mloc, schur_nloc, schur_lld;
    MUMPS_INT      mblock, nblock, nprow, npcol;
    MUMPS_INT      ld_rhsintr;
    MUMPS_INT      info[80],infog[80];
    XMUMPS_REAL    rinfo[40], rinfog[40];

    /* Null space */
    MUMPS_INT      deficiency;
    MUMPS_INT      *pivnul_list;
    MUMPS_INT      *mapping;
    XMUMPS_REAL    *singular_values;

    /* Schur */
    MUMPS_INT      size_schur;
    MUMPS_INT      *listvar_schur;
    XMUMPS_COMPLEX *schur;

    /* user workspace */
    XMUMPS_COMPLEX *wk_user;

    /* Version number: length=30 in FORTRAN + 1 for final \0 + 1 for alignment */
    char version_number[MUMPS_VERSION_MAX_LEN + 1 + 1];
    /* For out-of-core */
    char ooc_tmpdir[1024];
    char ooc_prefix[256];
    /* To save the matrix in matrix market format */
    char write_problem[1024];
    MUMPS_INT      lwk_user;
    /* For save/restore feature */
    char save_dir[1024];
    char save_prefix[256];

    /* Metis options */
    MUMPS_INT metis_options[40];

    /* Internal parameters */
    MUMPS_INT      instance_number;
} XMUMPS_STRUC_C;


void MUMPS_CALL
xmumps_c( XMUMPS_STRUC_C * xmumps_par );

#ifdef __cplusplus
}
#endif

#endif /* XMUMPS_C_H */
