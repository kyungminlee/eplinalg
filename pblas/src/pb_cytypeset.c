/* ---------------------------------------------------------------------
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  Include files
*/
#include "pblas.h"
#include "PBpblas.h"
#include "PBtools.h"
#include "PBblacs.h"
#include "PBblas.h"

PBTYP_T * PB_Cytypeset(void)
{
/*
*  Purpose
*  =======
*
*  PB_Cytypeset on the first call initializes a static structure contai-
*  ning typed information and returns a pointer to it.  The  other calls
*  to this routine just returns this pointer.
*
*  -- Written on April 1, 1998 by
*     R. Clint Whaley, University of Tennessee, Knoxville 37996, USA.
*
*  ---------------------------------------------------------------------
*/
/*
*  .. Local Scalars ..
*/
   static Int     setup=0;
   static PBTYP_T TypeStruct;
   static cmplxE zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = DCPLX;
   TypeStruct.usiz = sizeof( EREAL  );
   TypeStruct.size = sizeof( cmplxE );

   zero  [REAL_PART] = ZERO;
   zero  [IMAG_PART] = ZERO;
   one   [REAL_PART] =  ONE;
   one   [IMAG_PART] = ZERO;
   negone[REAL_PART] = -ONE;
   negone[IMAG_PART] = ZERO;

   TypeStruct.zero      = ((char *) zero);
   TypeStruct.one       = ((char *) one);
   TypeStruct.negone    = ((char *) negone);

   TypeStruct.Cgesd2d   = Cygesd2d;
   TypeStruct.Cgerv2d   = Cygerv2d;
   TypeStruct.Cgebs2d   = Cygebs2d;
   TypeStruct.Cgebr2d   = Cygebr2d;
   TypeStruct.Cgsum2d   = Cygsum2d;

   TypeStruct.Fmmadd    = ymmadd_;
   TypeStruct.Fmmcadd   = ymmcadd_;
   TypeStruct.Fmmtadd   = ymmtadd_;
   TypeStruct.Fmmtcadd  = ymmtcadd_;
   TypeStruct.Fmmdda    = ymmdda_;
   TypeStruct.Fmmddac   = ymmddac_;
   TypeStruct.Fmmddat   = ymmddat_;
   TypeStruct.Fmmddact  = ymmddact_;

   TypeStruct.Fcshft    = ycshft_;
   TypeStruct.Frshft    = yrshft_;

   TypeStruct.Fvvdotu   = yvvdotu_;
   TypeStruct.Fvvdotc   = yvvdotc_;

   TypeStruct.Fset      = yset_;

   TypeStruct.Ftzpad    = ytzpad_;
   TypeStruct.Ftzpadcpy = ytzpadcpy_;
   TypeStruct.Ftzscal   = ytzscal_;
   TypeStruct.Fhescal   = yhescal_;
   TypeStruct.Ftzcnjg   = ytzcnjg_;

   TypeStruct.Faxpy     = yaxpy_;
   TypeStruct.Fcopy     = ycopy_;
   TypeStruct.Fswap     = yswap_;

   TypeStruct.Fgemv     = ygemv_;
   TypeStruct.Fsymv     = ysymv_;
   TypeStruct.Fhemv     = yhemv_;
   TypeStruct.Ftrmv     = ytrmv_;
   TypeStruct.Ftrsv     = ytrsv_;
   TypeStruct.Fagemv    = yagemv_;
   TypeStruct.Fasymv    = yasymv_;
   TypeStruct.Fahemv    = yahemv_;
   TypeStruct.Fatrmv    = yatrmv_;

   TypeStruct.Fgerc     = ygerc_;
   TypeStruct.Fgeru     = ygeru_;
   TypeStruct.Fsyr      = ysyr_;
   TypeStruct.Fher      = yher_;
   TypeStruct.Fsyr2     = ysyr2_;
   TypeStruct.Fher2     = yher2_;

   TypeStruct.Fgemm     = ygemm_;
   TypeStruct.Fsymm     = ysymm_;
   TypeStruct.Fhemm     = yhemm_;
   TypeStruct.Fsyrk     = ysyrk_;
   TypeStruct.Fherk     = yherk_;
   TypeStruct.Fsyr2k    = ysyr2k_;
   TypeStruct.Fher2k    = yher2k_;
   TypeStruct.Ftrmm     = ytrmm_;
   TypeStruct.Ftrsm     = ytrsm_;

   return( &TypeStruct );
/*
*  End of PB_Cytypeset
*/
}
