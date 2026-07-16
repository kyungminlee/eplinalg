/*
 * Extended-precision numeric widths for the migrated MUMPS *bridge*
 * (x87 80-bit ``long double`` (REAL/COMPLEX(KIND=10))).
 *
 * Force-included (``-include mumps_c_types_extended.h``) ONLY into the
 * upstream ``mumps_c.c`` bridge objects — NOT into the standalone
 * emumps_c.h / ymumps_c.h consumer headers. It ``#include``s
 * upstream ``mumps_c_types.h`` (the bridge needs ``MUMPS_INT`` and the
 * ``MUMPS_ARITH_*`` selector flags), then overrides upstream's plain
 * ``double`` ``D/ZMUMPS_*`` widths so the upstream ``dmumps_c.h`` /
 * ``zmumps_c.h`` struct body — renamed into the bridge via
 * -DEMUMPS_STRUC_C… — lays out at the extended width, keeping the bridge
 * object ABI-identical to the consumer struct.
 *
 * Derived from MUMPS 5.8.2 (Copyright 1991-2026 CERFACS, CNRS, ENS Lyon,
 * INP Toulouse, Inria, Mumps Technologies, University of Bordeaux),
 * released under the CeCILL-C license.
 */

#ifndef MUMPS_C_TYPES_EXTENDED_H
#define MUMPS_C_TYPES_EXTENDED_H

/* Upstream MUMPS_INT / MUMPS_INT8 / MUMPS_ARITH_* / mumps_ftnlen and the
 * plain-``double`` S/C/D/ZMUMPS_* widths (the last of which we override
 * below). */
#include "mumps_c_types.h"

/* Extended real / complex widths as named types. */
#ifndef MUMPS_LONG_DOUBLE_TYPES
#define MUMPS_LONG_DOUBLE_TYPES
typedef long double                  mumps_long_double;
typedef struct { long double r, i; } mumps_long_double_complex;
#endif

/* Override upstream's plain-``double`` widths so the upstream C bridge
 * (mumps_c.c, compiled with -DDMUMPS_STRUC_C=EMUMPS_STRUC_C etc. and
 * force-including this header) instantiates the SAME layout the
 * emumps_c.h / ymumps_c.h headers declare. */
#undef  DMUMPS_REAL
#undef  DMUMPS_COMPLEX
#undef  ZMUMPS_REAL
#undef  ZMUMPS_COMPLEX
#define DMUMPS_REAL    mumps_long_double
#define DMUMPS_COMPLEX mumps_long_double
#define ZMUMPS_REAL    mumps_long_double
#define ZMUMPS_COMPLEX mumps_long_double_complex

/* --- Fortran-callback symbol migration ----------------------------
 * The scaling / RHS / singular-value callbacks are defined in mumps_c.c
 * (F_SYM_ARITH(assign_colsca,ASSIGN_COLSCA) -> with -DAdd_ the pasted
 * token dmumps_assign_colsca_, and the zmumps_ analogues). The migrated
 * *.F now CALL these at the target prefix (emumps_.../ymumps_...,
 * forced via recipes/mumps.yaml extra_renames), so retarget the bridge's
 * own definitions to match -- otherwise the extended build would leak
 * plain-double d/z callback symbols and collide with the genuine
 * dmumps_c/zmumps_c bridge. A ## paste is re-scanned, so a plain #define
 * on the final symbol rebinds it (same trick as -Ddmumps_set_tmp_ptr_=).
 * This header is force-included into BOTH bridge objects; the unused
 * arithmetic's set is simply never pasted, hence inert. */

#define dmumps_assign_colsca_                emumps_assign_colsca_
#define dmumps_assign_colsca_loc_            emumps_assign_colsca_loc_
#define dmumps_assign_rhsintr_               emumps_assign_rhsintr_
#define dmumps_assign_rowsca_                emumps_assign_rowsca_
#define dmumps_assign_rowsca_loc_            emumps_assign_rowsca_loc_
#define dmumps_assign_singular_values_       emumps_assign_singular_values_
#define dmumps_nullify_c_colsca_             emumps_nullify_c_colsca_
#define dmumps_nullify_c_colsca_loc_         emumps_nullify_c_colsca_loc_
#define dmumps_nullify_c_rhsintr_            emumps_nullify_c_rhsintr_
#define dmumps_nullify_c_rowsca_             emumps_nullify_c_rowsca_
#define dmumps_nullify_c_rowsca_loc_         emumps_nullify_c_rowsca_loc_
#define dmumps_nullify_c_sing_values_        emumps_nullify_c_sing_values_
#define dmumps_set_tmp_ptr_c_                emumps_set_tmp_ptr_c_
#define zmumps_assign_colsca_                ymumps_assign_colsca_
#define zmumps_assign_colsca_loc_            ymumps_assign_colsca_loc_
#define zmumps_assign_rhsintr_               ymumps_assign_rhsintr_
#define zmumps_assign_rowsca_                ymumps_assign_rowsca_
#define zmumps_assign_rowsca_loc_            ymumps_assign_rowsca_loc_
#define zmumps_assign_singular_values_       ymumps_assign_singular_values_
#define zmumps_nullify_c_colsca_             ymumps_nullify_c_colsca_
#define zmumps_nullify_c_colsca_loc_         ymumps_nullify_c_colsca_loc_
#define zmumps_nullify_c_rhsintr_            ymumps_nullify_c_rhsintr_
#define zmumps_nullify_c_rowsca_             ymumps_nullify_c_rowsca_
#define zmumps_nullify_c_rowsca_loc_         ymumps_nullify_c_rowsca_loc_
#define zmumps_nullify_c_sing_values_        ymumps_nullify_c_sing_values_
#define zmumps_set_tmp_ptr_c_                ymumps_set_tmp_ptr_c_

#endif /* MUMPS_C_TYPES_EXTENDED_H */
