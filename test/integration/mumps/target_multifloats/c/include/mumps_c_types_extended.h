/*
 * Extended-precision numeric widths for the migrated MUMPS *bridge*
 * (double-double multifloats (TYPE(real64x2)/TYPE(cmplx64x2))).
 *
 * Force-included (``-include mumps_c_types_extended.h``) ONLY into the
 * upstream ``mumps_c.c`` bridge objects — NOT into the standalone
 * mmumps_c.h / wmumps_c.h consumer headers. It ``#include``s
 * upstream ``mumps_c_types.h`` (the bridge needs ``MUMPS_INT`` and the
 * ``MUMPS_ARITH_*`` selector flags), then overrides upstream's plain
 * ``double`` ``D/ZMUMPS_*`` widths so the upstream ``dmumps_c.h`` /
 * ``zmumps_c.h`` struct body — renamed into the bridge via
 * -DMMUMPS_STRUC_C… — lays out at the extended width, keeping the bridge
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
#ifndef MUMPS_FLOAT64X2_TYPES
#define MUMPS_FLOAT64X2_TYPES
typedef struct { double limbs[2]; }      mumps_float64x2;
typedef struct { mumps_float64x2 r, i; } mumps_complex64x2;
#endif

/* Override upstream's plain-``double`` widths so the upstream C bridge
 * (mumps_c.c, compiled with -DDMUMPS_STRUC_C=MMUMPS_STRUC_C etc. and
 * force-including this header) instantiates the SAME layout the
 * mmumps_c.h / wmumps_c.h headers declare. */
#undef  DMUMPS_REAL
#undef  DMUMPS_COMPLEX
#undef  ZMUMPS_REAL
#undef  ZMUMPS_COMPLEX
#define DMUMPS_REAL    mumps_float64x2
#define DMUMPS_COMPLEX mumps_float64x2
#define ZMUMPS_REAL    mumps_float64x2
#define ZMUMPS_COMPLEX mumps_complex64x2

/* --- Fortran-callback symbol migration ----------------------------
 * The scaling / RHS / singular-value callbacks are defined in mumps_c.c
 * (F_SYM_ARITH(assign_colsca,ASSIGN_COLSCA) -> with -DAdd_ the pasted
 * token dmumps_assign_colsca_, and the zmumps_ analogues). The migrated
 * *.F now CALL these at the target prefix (mmumps_.../wmumps_...,
 * forced via codegen/recipes/mumps.yaml extra_renames), so retarget the bridge's
 * own definitions to match -- otherwise the extended build would leak
 * plain-double d/z callback symbols and collide with the genuine
 * dmumps_c/zmumps_c bridge. A ## paste is re-scanned, so a plain #define
 * on the final symbol rebinds it (same trick as -Ddmumps_set_tmp_ptr_=).
 * This header is force-included into BOTH bridge objects; the unused
 * arithmetic's set is simply never pasted, hence inert. */

#define dmumps_assign_colsca_                mmumps_assign_colsca_
#define dmumps_assign_colsca_loc_            mmumps_assign_colsca_loc_
#define dmumps_assign_rhsintr_               mmumps_assign_rhsintr_
#define dmumps_assign_rowsca_                mmumps_assign_rowsca_
#define dmumps_assign_rowsca_loc_            mmumps_assign_rowsca_loc_
#define dmumps_assign_singular_values_       mmumps_assign_singular_values_
#define dmumps_nullify_c_colsca_             mmumps_nullify_c_colsca_
#define dmumps_nullify_c_colsca_loc_         mmumps_nullify_c_colsca_loc_
#define dmumps_nullify_c_rhsintr_            mmumps_nullify_c_rhsintr_
#define dmumps_nullify_c_rowsca_             mmumps_nullify_c_rowsca_
#define dmumps_nullify_c_rowsca_loc_         mmumps_nullify_c_rowsca_loc_
#define dmumps_nullify_c_sing_values_        mmumps_nullify_c_sing_values_
#define dmumps_set_tmp_ptr_c_                mmumps_set_tmp_ptr_c_
#define zmumps_assign_colsca_                wmumps_assign_colsca_
#define zmumps_assign_colsca_loc_            wmumps_assign_colsca_loc_
#define zmumps_assign_rhsintr_               wmumps_assign_rhsintr_
#define zmumps_assign_rowsca_                wmumps_assign_rowsca_
#define zmumps_assign_rowsca_loc_            wmumps_assign_rowsca_loc_
#define zmumps_assign_singular_values_       wmumps_assign_singular_values_
#define zmumps_nullify_c_colsca_             wmumps_nullify_c_colsca_
#define zmumps_nullify_c_colsca_loc_         wmumps_nullify_c_colsca_loc_
#define zmumps_nullify_c_rhsintr_            wmumps_nullify_c_rhsintr_
#define zmumps_nullify_c_rowsca_             wmumps_nullify_c_rowsca_
#define zmumps_nullify_c_rowsca_loc_         wmumps_nullify_c_rowsca_loc_
#define zmumps_nullify_c_sing_values_        wmumps_nullify_c_sing_values_
#define zmumps_set_tmp_ptr_c_                wmumps_set_tmp_ptr_c_

#endif /* MUMPS_C_TYPES_EXTENDED_H */
