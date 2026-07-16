/* pb_cmtypeset.c -- multifloats replacement for PB_Cdtypeset.
 *
 * The migrated-by-regex version cannot initialize the static
 * zero/one/negone constants with C operators on the float64x2
 * struct (e.g. ``zero = ZERO`` where ZERO is the macro 0.0). This
 * hand-written version uses compound literals from libmfc to set
 * the constants and points the function pointers at the migrated
 * BLACS / BLAS / PBBLAS / PTZBLAS routines that the multifloats
 * BLACS migration produced.
 */
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"
#include "multifloats_bridge.h"


/*
*  Hidden-CHARACTER-length ABI trampolines.  Each F<slot> callback
*  below is reached by the type-generic drivers through a C prototype
*  that passes NO Fortran CHARACTER lengths; the migrated gfortran
*  leaf expects them and spills scratch into the (unprovisioned)
*  hidden-length slots, clobbering the caller frame.  These
*  trampolines carry the driver's no-length signature inward and
*  append the lengths (all 1) outward.  See PBcharshim.h.
*/
#include "PBcharshim.h"

EP_MK_TZPAD(tzpad, mtzpad_)
EP_MK_TZPADCPY(tzpadcpy, mtzpadcpy_)
EP_MK_TZSCAL(tzscal, mtzscal_)
EP_MK_TZSCAL(hescal, mtzscal_)
EP_MK_TZSCAL(tzcnjg, mtzscal_)
EP_MK_GEMV(gemv, mgemv_)
EP_MK_SYMV(symv, msymv_)
EP_MK_SYMV(hemv, msymv_)
EP_MK_TRMV(trmv, mtrmv_)
EP_MK_TRMV(trsv, mtrsv_)
EP_MK_GEMV(agemv, magemv_)
EP_MK_SYMV(asymv, masymv_)
EP_MK_SYMV(ahemv, masymv_)
EP_MK_ATRMV(atrmv, matrmv_)
EP_MK_SYR(syr, msyr_)
EP_MK_SYR(her, msyr_)
EP_MK_SYR2(syr2, msyr2_)
EP_MK_SYR2(her2, msyr2_)
EP_MK_GEMM(gemm, mgemm_)
EP_MK_SYMM(symm, msymm_)
EP_MK_SYMM(hemm, msymm_)
EP_MK_SYRK(syrk, msyrk_)
EP_MK_SYRK(herk, msyrk_)
EP_MK_SYR2K(syr2k, msyr2k_)
EP_MK_SYR2K(her2k, msyr2k_)
EP_MK_TRMM(trmm, mtrmm_)
EP_MK_TRMM(trsm, mtrsm_)

PBTYP_T * PB_Cmtypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static float64x2 zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DREAL;
   TypeStruct.usiz = sizeof(float64x2);
   TypeStruct.size = sizeof(float64x2);

   zero   = 0.0;
   one    = 1.0;
   negone = -1.0;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cmgesd2d;
   TypeStruct.Cgerv2d   = Cmgerv2d;
   TypeStruct.Cgebs2d   = Cmgebs2d;
   TypeStruct.Cgebr2d   = Cmgebr2d;
   TypeStruct.Cgsum2d   = Cmgsum2d;

   TypeStruct.Fmmadd    = mmmadd_;
   TypeStruct.Fmmcadd   = mmmcadd_;
   TypeStruct.Fmmtadd   = mmmtadd_;
   TypeStruct.Fmmtcadd  = mmmtcadd_;
   TypeStruct.Fmmdda    = mmmdda_;
   TypeStruct.Fmmddac   = mmmddac_;
   TypeStruct.Fmmddat   = mmmddat_;
   TypeStruct.Fmmddact  = mmmddact_;

   TypeStruct.Fcshft    = mcshft_;
   TypeStruct.Frshft    = mrshft_;

   TypeStruct.Fvvdotu   = mvvdot_;
   TypeStruct.Fvvdotc   = mvvdot_;

   TypeStruct.Fset      = mset_;

   TypeStruct.Ftzpad    = (TZPAD_T) EP_SHIM(tzpad);
   TypeStruct.Ftzpadcpy = (TZPADCPY_T) EP_SHIM(tzpadcpy);
   TypeStruct.Ftzscal   = (TZSCAL_T) EP_SHIM(tzscal);
   TypeStruct.Fhescal   = (TZSCAL_T) EP_SHIM(hescal);
   TypeStruct.Ftzcnjg   = (TZSCAL_T) EP_SHIM(tzcnjg);

   TypeStruct.Faxpy     = maxpy_;
   TypeStruct.Fcopy     = mcopy_;
   TypeStruct.Fswap     = mswap_;

   TypeStruct.Fgemv     = (GEMV_T) EP_SHIM(gemv);
   TypeStruct.Fsymv     = (SYMV_T) EP_SHIM(symv);
   TypeStruct.Fhemv     = (HEMV_T) EP_SHIM(hemv);
   TypeStruct.Ftrmv     = (TRMV_T) EP_SHIM(trmv);
   TypeStruct.Ftrsv     = (TRSV_T) EP_SHIM(trsv);
   TypeStruct.Fagemv    = (AGEMV_T) EP_SHIM(agemv);
   TypeStruct.Fasymv    = (ASYMV_T) EP_SHIM(asymv);
   TypeStruct.Fahemv    = (AHEMV_T) EP_SHIM(ahemv);
   TypeStruct.Fatrmv    = (ATRMV_T) EP_SHIM(atrmv);

   TypeStruct.Fgerc     = mger_;
   TypeStruct.Fgeru     = mger_;
   TypeStruct.Fsyr      = (SYR_T) EP_SHIM(syr);
   TypeStruct.Fher      = (HER_T) EP_SHIM(her);
   TypeStruct.Fsyr2     = (SYR2_T) EP_SHIM(syr2);
   TypeStruct.Fher2     = (HER2_T) EP_SHIM(her2);

   TypeStruct.Fgemm     = (GEMM_T) EP_SHIM(gemm);
   TypeStruct.Fsymm     = (SYMM_T) EP_SHIM(symm);
   TypeStruct.Fhemm     = (HEMM_T) EP_SHIM(hemm);
   TypeStruct.Fsyrk     = (SYRK_T) EP_SHIM(syrk);
   TypeStruct.Fherk     = (HERK_T) EP_SHIM(herk);
   TypeStruct.Fsyr2k    = (SYR2K_T) EP_SHIM(syr2k);
   TypeStruct.Fher2k    = (HER2K_T) EP_SHIM(her2k);
   TypeStruct.Ftrmm     = (TRMM_T) EP_SHIM(trmm);
   TypeStruct.Ftrsm     = (TRSM_T) EP_SHIM(trsm);

   return( &TypeStruct );
}
