/* pb_cwtypeset.c -- multifloats replacement for PB_Cztypeset.
 *
 * Mirror of pb_cmtypeset.c for the complex double-double variant.
 * Initializes the static zero/one/negone constants in array form
 * (cmplxDD = float64x2[2]) and points the function pointers at
 * the migrated complex BLACS / BLAS / PBBLAS / PTZBLAS routines.
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

EP_MK_TZPAD(tzpad, wtzpad_)
EP_MK_TZPADCPY(tzpadcpy, wtzpadcpy_)
EP_MK_TZSCAL(tzscal, wtzscal_)
EP_MK_TZSCAL(hescal, whescal_)
EP_MK_TZSCAL(tzcnjg, wtzcnjg_)
EP_MK_GEMV(gemv, wgemv_)
EP_MK_SYMV(symv, wsymv_)
EP_MK_SYMV(hemv, whemv_)
EP_MK_TRMV(trmv, wtrmv_)
EP_MK_TRMV(trsv, wtrsv_)
EP_MK_GEMV(agemv, wagemv_)
EP_MK_SYMV(asymv, wasymv_)
EP_MK_SYMV(ahemv, wahemv_)
EP_MK_ATRMV(atrmv, watrmv_)
EP_MK_SYR(syr, wsyr_)
EP_MK_SYR(her, wher_)
EP_MK_SYR2(syr2, wsyr2_)
EP_MK_SYR2(her2, wher2_)
EP_MK_GEMM(gemm, wgemm_)
EP_MK_SYMM(symm, wsymm_)
EP_MK_SYMM(hemm, whemm_)
EP_MK_SYRK(syrk, wsyrk_)
EP_MK_SYRK(herk, wherk_)
EP_MK_SYR2K(syr2k, wsyr2k_)
EP_MK_SYR2K(her2k, wher2k_)
EP_MK_TRMM(trmm, wtrmm_)
EP_MK_TRMM(trsm, wtrsm_)

PBTYP_T * PB_Cwtypeset(void)
{
   static Int     setup = 0;
   static PBTYP_T TypeStruct;
   static cmplxDD zero, one, negone;

   if( setup ) return( &TypeStruct );
   setup = 1;

   TypeStruct.type = DCPLX;
   TypeStruct.usiz = sizeof(float64x2);
   TypeStruct.size = sizeof(cmplxDD);

   zero[REAL_PART]   = 0.0;
   zero[IMAG_PART]   = 0.0;
   one[REAL_PART]    = 1.0;
   one[IMAG_PART]    = 0.0;
   negone[REAL_PART] = -1.0;
   negone[IMAG_PART] = 0.0;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Cwgesd2d;
   TypeStruct.Cgerv2d   = Cwgerv2d;
   TypeStruct.Cgebs2d   = Cwgebs2d;
   TypeStruct.Cgebr2d   = Cwgebr2d;
   TypeStruct.Cgsum2d   = Cwgsum2d;

   TypeStruct.Fmmadd    = wmmadd_;
   TypeStruct.Fmmcadd   = wmmcadd_;
   TypeStruct.Fmmtadd   = wmmtadd_;
   TypeStruct.Fmmtcadd  = wmmtcadd_;
   TypeStruct.Fmmdda    = wmmdda_;
   TypeStruct.Fmmddac   = wmmddac_;
   TypeStruct.Fmmddat   = wmmddat_;
   TypeStruct.Fmmddact  = wmmddact_;

   TypeStruct.Fcshft    = wcshft_;
   TypeStruct.Frshft    = wrshft_;

   TypeStruct.Fvvdotu   = wvvdotu_;
   TypeStruct.Fvvdotc   = wvvdotc_;

   TypeStruct.Fset      = wset_;

   TypeStruct.Ftzpad    = (TZPAD_T) EP_SHIM(tzpad);
   TypeStruct.Ftzpadcpy = (TZPADCPY_T) EP_SHIM(tzpadcpy);
   TypeStruct.Ftzscal   = (TZSCAL_T) EP_SHIM(tzscal);
   TypeStruct.Fhescal   = (TZSCAL_T) EP_SHIM(hescal);
   TypeStruct.Ftzcnjg   = (TZSCAL_T) EP_SHIM(tzcnjg);

   TypeStruct.Faxpy     = waxpy_;
   TypeStruct.Fcopy     = wcopy_;
   TypeStruct.Fswap     = wswap_;

   TypeStruct.Fgemv     = (GEMV_T) EP_SHIM(gemv);
   TypeStruct.Fsymv     = (SYMV_T) EP_SHIM(symv);
   TypeStruct.Fhemv     = (HEMV_T) EP_SHIM(hemv);
   TypeStruct.Ftrmv     = (TRMV_T) EP_SHIM(trmv);
   TypeStruct.Ftrsv     = (TRSV_T) EP_SHIM(trsv);
   TypeStruct.Fagemv    = (AGEMV_T) EP_SHIM(agemv);
   TypeStruct.Fasymv    = (ASYMV_T) EP_SHIM(asymv);
   TypeStruct.Fahemv    = (AHEMV_T) EP_SHIM(ahemv);
   TypeStruct.Fatrmv    = (ATRMV_T) EP_SHIM(atrmv);

   TypeStruct.Fgerc     = wgerc_;
   TypeStruct.Fgeru     = wgeru_;
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
