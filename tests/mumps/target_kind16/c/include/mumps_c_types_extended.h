/*
 * Extended-precision numeric widths for the migrated MUMPS *bridge*
 * (``__float128`` (REAL/COMPLEX(KIND=16))).
 *
 * Force-included (``-include mumps_c_types_extended.h``) ONLY into the
 * upstream ``mumps_c.c`` bridge objects — NOT into the standalone
 * qmumps_c.h / xmumps_c.h consumer headers. It ``#include``s
 * upstream ``mumps_c_types.h`` (the bridge needs ``MUMPS_INT`` and the
 * ``MUMPS_ARITH_*`` selector flags), then overrides upstream's plain
 * ``double`` ``D/ZMUMPS_*`` widths so the upstream ``dmumps_c.h`` /
 * ``zmumps_c.h`` struct body — renamed into the bridge via
 * -DQMUMPS_STRUC_C… — lays out at the extended width, keeping the bridge
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
#ifndef MUMPS_FLOAT128_TYPES
#define MUMPS_FLOAT128_TYPES
typedef __float128                  mumps_float128;
typedef struct { __float128 r, i; } mumps_float128_complex;
#endif

/* Override upstream's plain-``double`` widths so the upstream C bridge
 * (mumps_c.c, compiled with -DDMUMPS_STRUC_C=QMUMPS_STRUC_C etc. and
 * force-including this header) instantiates the SAME layout the
 * qmumps_c.h / xmumps_c.h headers declare. */
#undef  DMUMPS_REAL
#undef  DMUMPS_COMPLEX
#undef  ZMUMPS_REAL
#undef  ZMUMPS_COMPLEX
#define DMUMPS_REAL    mumps_float128
#define DMUMPS_COMPLEX mumps_float128
#define ZMUMPS_REAL    mumps_float128
#define ZMUMPS_COMPLEX mumps_float128_complex

/* --- Fortran-callback symbol migration ----------------------------
 * The scaling / RHS / singular-value callbacks are defined in mumps_c.c
 * (F_SYM_ARITH(assign_colsca,ASSIGN_COLSCA) -> with -DAdd_ the pasted
 * token dmumps_assign_colsca_, and the zmumps_ analogues). The migrated
 * *.F now CALL these at the target prefix (qmumps_.../xmumps_...,
 * forced via recipes/mumps.yaml extra_renames), so retarget the bridge's
 * own definitions to match -- otherwise the extended build would leak
 * plain-double d/z callback symbols and collide with the genuine
 * dmumps_c/zmumps_c bridge. A ## paste is re-scanned, so a plain #define
 * on the final symbol rebinds it (same trick as -Ddmumps_set_tmp_ptr_=).
 * This header is force-included into BOTH bridge objects; the unused
 * arithmetic's set is simply never pasted, hence inert. */

#define dmumps_assign_colsca_                qmumps_assign_colsca_
#define dmumps_assign_colsca_loc_            qmumps_assign_colsca_loc_
#define dmumps_assign_rhsintr_               qmumps_assign_rhsintr_
#define dmumps_assign_rowsca_                qmumps_assign_rowsca_
#define dmumps_assign_rowsca_loc_            qmumps_assign_rowsca_loc_
#define dmumps_assign_singular_values_       qmumps_assign_singular_values_
#define dmumps_nullify_c_colsca_             qmumps_nullify_c_colsca_
#define dmumps_nullify_c_colsca_loc_         qmumps_nullify_c_colsca_loc_
#define dmumps_nullify_c_rhsintr_            qmumps_nullify_c_rhsintr_
#define dmumps_nullify_c_rowsca_             qmumps_nullify_c_rowsca_
#define dmumps_nullify_c_rowsca_loc_         qmumps_nullify_c_rowsca_loc_
#define dmumps_nullify_c_sing_values_        qmumps_nullify_c_sing_values_
#define dmumps_set_tmp_ptr_c_                qmumps_set_tmp_ptr_c_
#define zmumps_assign_colsca_                xmumps_assign_colsca_
#define zmumps_assign_colsca_loc_            xmumps_assign_colsca_loc_
#define zmumps_assign_rhsintr_               xmumps_assign_rhsintr_
#define zmumps_assign_rowsca_                xmumps_assign_rowsca_
#define zmumps_assign_rowsca_loc_            xmumps_assign_rowsca_loc_
#define zmumps_assign_singular_values_       xmumps_assign_singular_values_
#define zmumps_nullify_c_colsca_             xmumps_nullify_c_colsca_
#define zmumps_nullify_c_colsca_loc_         xmumps_nullify_c_colsca_loc_
#define zmumps_nullify_c_rhsintr_            xmumps_nullify_c_rhsintr_
#define zmumps_nullify_c_rowsca_             xmumps_nullify_c_rowsca_
#define zmumps_nullify_c_rowsca_loc_         xmumps_nullify_c_rowsca_loc_
#define zmumps_nullify_c_sing_values_        xmumps_nullify_c_sing_values_
#define zmumps_set_tmp_ptr_c_                xmumps_set_tmp_ptr_c_

#endif /* MUMPS_C_TYPES_EXTENDED_H */
