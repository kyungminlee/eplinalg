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

PBTYP_T * PB_Cetypeset(void)
{
/*
*  Purpose
*  =======
*
*  PB_Cetypeset on the first call initializes a static structure contai-
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
   static EREAL  zero, one, negone;
/* ..
*  .. Executable Statements ..
*
*/
   if( setup ) return( &TypeStruct );

   setup = 1;

   TypeStruct.type = DREAL;
   TypeStruct.usiz = sizeof( EREAL );
   TypeStruct.size = sizeof( EREAL );

   zero   = ZERO;
   one    =  ONE;
   negone = -ONE;

   TypeStruct.zero      = (char *) (&zero);
   TypeStruct.one       = (char *) (&one);
   TypeStruct.negone    = (char *) (&negone);

   TypeStruct.Cgesd2d   = Cegesd2d;
   TypeStruct.Cgerv2d   = Cegerv2d;
   TypeStruct.Cgebs2d   = Cegebs2d;
   TypeStruct.Cgebr2d   = Cegebr2d;
   TypeStruct.Cgsum2d   = Cegsum2d;

   TypeStruct.Fmmadd    = emmadd_;
   TypeStruct.Fmmcadd   = emmcadd_;
   TypeStruct.Fmmtadd   = emmtadd_;
   TypeStruct.Fmmtcadd  = emmtcadd_;
   TypeStruct.Fmmdda    = emmdda_;
   TypeStruct.Fmmddac   = emmddac_;
   TypeStruct.Fmmddat   = emmddat_;
   TypeStruct.Fmmddact  = emmddact_;

   TypeStruct.Fcshft    = ecshft_;
   TypeStruct.Frshft    = ershft_;

   TypeStruct.Fvvdotu   = evvdot_;
   TypeStruct.Fvvdotc   = evvdot_;

   TypeStruct.Fset      = eset_;

   TypeStruct.Ftzpad    = etzpad_;
   TypeStruct.Ftzpadcpy = etzpadcpy_;
   TypeStruct.Ftzscal   = etzscal_;
   TypeStruct.Fhescal   = etzscal_;
   TypeStruct.Ftzcnjg   = etzscal_;

   TypeStruct.Faxpy     = eaxpy_;
   TypeStruct.Fcopy     = ecopy_;
   TypeStruct.Fswap     = eswap_;

   TypeStruct.Fgemv     = egemv_;
   TypeStruct.Fsymv     = esymv_;
   TypeStruct.Fhemv     = esymv_;
   TypeStruct.Ftrmv     = etrmv_;
   TypeStruct.Ftrsv     = etrsv_;
   TypeStruct.Fagemv    = eagemv_;
   TypeStruct.Fasymv    = easymv_;
   TypeStruct.Fahemv    = easymv_;
   TypeStruct.Fatrmv    = eatrmv_;

   TypeStruct.Fgerc     = eger_;
   TypeStruct.Fgeru     = eger_;
   TypeStruct.Fsyr      = esyr_;
   TypeStruct.Fher      = esyr_;
   TypeStruct.Fsyr2     = esyr2_;
   TypeStruct.Fher2     = esyr2_;

   TypeStruct.Fgemm     = egemm_;
   TypeStruct.Fsymm     = esymm_;
   TypeStruct.Fhemm     = esymm_;
   TypeStruct.Fsyrk     = esyrk_;
   TypeStruct.Fherk     = esyrk_;
   TypeStruct.Fsyr2k    = esyr2k_;
   TypeStruct.Fher2k    = esyr2k_;
   TypeStruct.Ftrmm     = etrmm_;
   TypeStruct.Ftrsm     = etrsm_;

   return( &TypeStruct );
/*
*  End of PB_Cetypeset
*/
}
