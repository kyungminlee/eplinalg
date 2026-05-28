/* ---------------------------------------------------------------------
*
*  -- ScaLAPACK routine (version 1.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     November 17, 1996
*
*  ---------------------------------------------------------------------
*/

/*
* This file includes the standard C libraries, as well as system
* dependent include files.  All PBLAS routines include this file.
*/

/*
* ========================================================================
* Machine Specific PBLAS macros
* ========================================================================
*/

#define _HAL_           0
#define _T3D_           1

#ifdef T3D
#define _MACH_          _T3D_
#endif

#ifndef _MACH_
#define _MACH_          _HAL_
#endif

/*
* ========================================================================
* Include files
* ========================================================================
*/
#include <stdio.h>
#include <stdlib.h>

#if( _MACH_ == _T3D_ )
#include <fortran.h>
#endif

/*
* ========================================================================
* FORTRAN <-> C interface
* ========================================================================
*
* These macros define how the PBLAS will be called. _F2C_ADD_ assumes 
* that they will be called by FORTRAN, which expects C routines to have
* an underscore postfixed to the name (Suns, and Intel machines expect
* this). _F2C_NOCHANGE indicates that FORTRAN will be calling, and that
* it expects the name called by FORTRAN to be identical to that compiled
* by the C (RS6K's do this).  _F2C_UPCASE says it expects C routines
* called by FORTRAN to be in all upcase (CRAY wants this).   
* _F2C_F77ISF2C indicates that the fortran "compiler" in use is
* actually f2c, a FORTRAN to C converter.
*/

#define _F2C_ADD_       0
#define _F2C_NOCHANGE   1
#define _F2C_UPCASE     2
#define _F2C_F77ISF2C   3

#ifdef UpCase
#define _F2C_CALL_      _F2C_UPCASE
#endif

#ifdef NoChange
#define _F2C_CALL_      _F2C_NOCHANGE
#endif

#ifdef Add_
#define _F2C_CALL_      _F2C_ADD_
#endif

#ifdef f77IsF2C
#define _F2C_CALL_      _F2C_F77ISF2C
#endif

#ifndef _F2C_CALL_
#define _F2C_CALL_      _F2C_ADD_
#endif

/*
* ========================================================================
* TYPE DEFINITIONS AND CONVERSION UTILITIES
* ========================================================================
*/

#ifndef Int
#define Int int
#endif

typedef long double EREAL;
typedef struct { EREAL re, im; } ECOMPLEX;

typedef struct { float  re, im; } complex;
typedef struct { double re, im; } complex16;

#if( _MACH_ == _T3D_ )
#define float  double
                       /* Type of character argument in a FORTRAN call */
#define F_CHAR          _fcd
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     ( _fcdtocp( (a) ) )
#define C2F_CHAR(a)     ( _cptofcd( (a), 1 ) )
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void   fortran                   /* Subroutine */
#define F_INTG_FCT      Int    fortran             /* INTEGER function */
#define F_DBLE_FCT      double fortran    /* DOUBLE PRECISION function */

#else
                       /* Type of character argument in a FORTRAN call */
typedef char *          F_CHAR;
                                     /* Character conversion utilities */
#define F2C_CHAR(a)     (a)
#define C2F_CHAR(a)     (a)
                                          /* Type of FORTRAN functions */
#define F_VOID_FCT      void                             /* Subroutine */
#define F_INTG_FCT      Int                        /* INTEGER function */
#define F_DBLE_FCT      double            /* DOUBLE PRECISION function */

#endif

/*
* ========================================================================
* #DEFINE MACRO CONSTANTS
* ========================================================================
*/
#define    DLEN_        9                     /* Length of a descriptor */
#define    DT_          0                     /* Descriptor Type        */
#define    CTXT_        1                              /* BLACS context */
#define    M_           2                      /* Global Number of Rows */
#define    N_           3                   /* Global Number of Columns */
#define    MB_          4                          /* Row Blocking Size */
#define    NB_          5                       /* Column Blocking Size */
#define    RSRC_        6                     /* Starting Processor Row */
#define    CSRC_        7                  /* Starting Processor Column */
#define    LLD_         8                    /* Local Leading Dimension */

/*
 * Descriptor types
 */
#define    BLOCK_CYCLIC_2D                1
#define    BLOCK_CYCLIC_INB_2D            2

#define    BROADCAST    "B"              /* Blacs operation definitions */
#define    COMBINE      "C"

#define    ALL          "A"                        /* Scope definitions */
#define    COLUMN       "C"
#define    ROW          "R"

#define    TOPDEF       " " /* Default BLACS topology, PB-BLAS routines */
#define    CTOPDEF      ' '
#define    TOPGET       "!"

#define    YES          "Y"
#define    NO           "N"

#define    MULLENFAC    2

#define    ONE          1.0
#define    ZERO         0.0

/*
* ========================================================================
* PREPROCESSOR MACRO FUNCTIONS USED FOR OPTIMIZATION & CONVENIENCE
* ========================================================================
*/

#define ABS(a)   (((a) < 0) ? -(a) : (a))

#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

#define CEIL(a,b) ( ((a)+(b)-1) / (b) )

#define Mlowcase(C) ( ((C) > 64 && (C) < 91) ? (C) | 32 : (C) )

#define Mupcase(C) ( ((C) > 96 && (C) < 123) ? (C) & 0xDF : (C) )

#define INDXG2L( iglob, nb, iproc, isrcproc, nprocs )\
    ( (nb) * ( ( (iglob)-1) / ( (nb) * (nprocs) ) ) +\
      ( ( (iglob) - 1 ) % (nb) ) + 1 )

#define INDXL2G( iloc, nb, iproc, isrcproc, nprocs )\
    ( (nprocs) * (nb) * ( ( (iloc) - 1 ) / (nb) ) +\
      ( ( (iloc) - 1 ) % (nb) ) +\
      ( ( (nprocs) + (iproc) - (isrcproc) ) % (nprocs) ) * (nb) + 1 )

#define INDXG2P( iglob, nb, iproc, isrcproc, nprocs ) \
    ( ( (isrcproc) + ( (iglob) - 1 ) / (nb) ) % (nprocs) )

#define MYROC0( nblocks, n, nb, nprocs )\
  ( ( (nblocks) % (nprocs) ) ? ( ( (nblocks) / (nprocs) ) * (nb) + (nb) )\
                   : ( ( (nblocks) / (nprocs) )* (nb) + ( (n) % (nb) ) ) )

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in).
* No redefinition necessary to have following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm_(...)
*           call pegemm(...)           void pegemm_(...)
*
* This is the default.
*/

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void PDGEMM(...)
*           call pegemm(...)           void PEGEMM(...)
*/
                                                            /* TOOLS */
#define ilcm_             ILCM
#define infog2l_          INFOG2L
#define numroc_           NUMROC
#define pstreecomb_       PSTREECOMB
#define petreecomb_       PETREECOMB
#define pdtreecomb_       PDTREECOMB
#define petreecomb_       PETREECOMB
#define pctreecomb_       PCTREECOMB
#define pytreecomb_       PYTREECOMB
#define pztreecomb_       PZTREECOMB
#define pytreecomb_       PYTREECOMB
#define scombamax_        SCOMBAMAX
#define ecombamax_        ECOMBAMAX
#define dcombamax_        DCOMBAMAX
#define ecombamax_        ECOMBAMAX
#define ccombamax_        CCOMBAMAX
#define ycombamax_        YCOMBAMAX
#define zcombamax_        ZCOMBAMAX
#define ycombamax_        YCOMBAMAX
#define scombnrm2_        SCOMBNRM2
#define ecombnrm2_        ECOMBNRM2
#define dcombnrm2_        DCOMBNRM2
#define ecombnrm2_        ECOMBNRM2

#define dlamov_           DLAMOV
#define elamov_           ELAMOV
#define slamov_           SLAMOV
#define elamov_           ELAMOV
#define clamov_           CLAMOV
#define ylamov_           YLAMOV
#define zlamov_           ZLAMOV
#define ylamov_           YLAMOV
#define dlacpy_           DLACPY
#define elacpy_           ELACPY
#define slacpy_           SLACPY
#define elacpy_           ELACPY
#define clacpy_           CLACPY
#define ylacpy_           YLACPY
#define zlacpy_           ZLACPY
#define ylacpy_           YLACPY
#define xerbla_           XERBLA
                                                            /* BLACS */
#define blacs_abort_      BLACS_ABORT
#define blacs_gridinfo_   BLACS_GRIDINFO

#define igesd2d_          IGESD2D
#define igebs2d_          IGEBS2D
#define itrsd2d_          ITRSD2D
#define itrbs2d_          ITRBS2D
#define igerv2d_          IGERV2D
#define igebr2d_          IGEBR2D
#define itrrv2d_          ITRRV2D
#define itrbr2d_          ITRBR2D
#define igamx2d_          IGAMX2D
#define igamn2d_          IGAMN2D
#define igsum2d_          IGSUM2D

#define sgesd2d_          SGESD2D
#define egesd2d_          EGESD2D
#define sgebs2d_          SGEBS2D
#define egebs2d_          EGEBS2D
#define strsd2d_          STRSD2D
#define etrsd2d_          ETRSD2D
#define strbs2d_          STRBS2D
#define etrbs2d_          ETRBS2D
#define sgerv2d_          SGERV2D
#define egerv2d_          EGERV2D
#define sgebr2d_          SGEBR2D
#define egebr2d_          EGEBR2D
#define strrv2d_          STRRV2D
#define etrrv2d_          ETRRV2D
#define strbr2d_          STRBR2D
#define etrbr2d_          ETRBR2D
#define sgamx2d_          SGAMX2D
#define egamx2d_          EGAMX2D
#define sgamn2d_          SGAMN2D
#define egamn2d_          EGAMN2D
#define sgsum2d_          SGSUM2D
#define egsum2d_          EGSUM2D

#define dgesd2d_          DGESD2D
#define egesd2d_          EGESD2D
#define dgebs2d_          DGEBS2D
#define egebs2d_          EGEBS2D
#define dtrsd2d_          DTRSD2D
#define etrsd2d_          ETRSD2D
#define dtrbs2d_          DTRBS2D
#define etrbs2d_          ETRBS2D
#define dgerv2d_          DGERV2D
#define egerv2d_          EGERV2D
#define dgebr2d_          DGEBR2D
#define egebr2d_          EGEBR2D
#define dtrrv2d_          DTRRV2D
#define etrrv2d_          ETRRV2D
#define dtrbr2d_          DTRBR2D
#define etrbr2d_          ETRBR2D
#define dgamx2d_          DGAMX2D
#define egamx2d_          EGAMX2D
#define dgamn2d_          DGAMN2D
#define egamn2d_          EGAMN2D
#define dgsum2d_          DGSUM2D
#define egsum2d_          EGSUM2D

#define cgesd2d_          CGESD2D
#define ygesd2d_          YGESD2D
#define cgebs2d_          CGEBS2D
#define ygebs2d_          YGEBS2D
#define ctrsd2d_          CTRSD2D
#define ytrsd2d_          YTRSD2D
#define ctrbs2d_          CTRBS2D
#define ytrbs2d_          YTRBS2D
#define cgerv2d_          CGERV2D
#define ygerv2d_          YGERV2D
#define cgebr2d_          CGEBR2D
#define ygebr2d_          YGEBR2D
#define ctrrv2d_          CTRRV2D
#define ytrrv2d_          YTRRV2D
#define ctrbr2d_          CTRBR2D
#define ytrbr2d_          YTRBR2D
#define cgamx2d_          CGAMX2D
#define ygamx2d_          YGAMX2D
#define cgamn2d_          CGAMN2D
#define ygamn2d_          YGAMN2D
#define cgsum2d_          CGSUM2D
#define ygsum2d_          YGSUM2D

#define zgesd2d_          ZGESD2D
#define ygesd2d_          YGESD2D
#define zgebs2d_          ZGEBS2D
#define ygebs2d_          YGEBS2D
#define ztrsd2d_          ZTRSD2D
#define ytrsd2d_          YTRSD2D
#define ztrbs2d_          ZTRBS2D
#define ytrbs2d_          YTRBS2D
#define zgerv2d_          ZGERV2D
#define ygerv2d_          YGERV2D
#define zgebr2d_          ZGEBR2D
#define ygebr2d_          YGEBR2D
#define ztrrv2d_          ZTRRV2D
#define ytrrv2d_          YTRRV2D
#define ztrbr2d_          ZTRBR2D
#define ytrbr2d_          YTRBR2D
#define zgamx2d_          ZGAMX2D
#define ygamx2d_          YGAMX2D
#define zgamn2d_          ZGAMN2D
#define ygamn2d_          YGAMN2D
#define zgsum2d_          ZGSUM2D
#define ygsum2d_          YGSUM2D
                                                     /* Level-1 BLAS */
#define srotg_            SROTG
#define erotg_            EROTG
#define srotmg_           SROTMG
#define erotmg_           EROTMG
#define srot_             SROT
#define erot_             EROT
#define srotm_            SROTM
#define erotm_            EROTM
#define sswap_            SSWAP
#define eswap_            ESWAP
#define sscal_            SSCAL
#define escal_            ESCAL
#define scopy_            SCOPY
#define ecopy_            ECOPY
#define saxpy_            SAXPY
#define eaxpy_            EAXPY
#define ssdot_            SSDOT
#define esdot_            ESDOT
#define isamax_           ISAMAX
#define ieamax_           IEAMAX

#define drotg_            DROTG
#define erotg_            EROTG
#define drotmg_           DROTMG
#define erotmg_           EROTMG
#define drot_             DROT
#define erot_             EROT
#define drotm_            DROTM
#define erotm_            EROTM
#define dswap_            DSWAP
#define eswap_            ESWAP
#define dscal_            DSCAL
#define escal_            ESCAL
#define dcopy_            DCOPY
#define ecopy_            ECOPY
#define daxpy_            DAXPY
#define eaxpy_            EAXPY
#define dddot_            DDDOT
#define dnrm2_            DNRM2
#define enrm2_            ENRM2
#define dsnrm2_           DSNRM2
#define dasum_            DASUM
#define easum_            EASUM
#define dsasum_           DSASUM
#define idamax_           IDAMAX
#define ieamax_           IEAMAX

#define cswap_            CSWAP
#define yswap_            YSWAP
#define cscal_            CSCAL
#define yscal_            YSCAL
#define csscal_           CSSCAL
#define yescal_           YESCAL
#define ccopy_            CCOPY
#define ycopy_            YCOPY
#define caxpy_            CAXPY
#define yaxpy_            YAXPY
#define ccdotu_           CCDOTU
#define yydotu_           YYDOTU
#define ccdotc_           CCDOTC
#define yydotc_           YYDOTC
#define icamax_           ICAMAX
#define iyamax_           IYAMAX

#define zswap_            ZSWAP
#define yswap_            YSWAP
#define zscal_            ZSCAL
#define yscal_            YSCAL
#define zdscal_           ZDSCAL
#define yescal_           YESCAL
#define zcopy_            ZCOPY
#define ycopy_            YCOPY
#define zaxpy_            ZAXPY
#define yaxpy_            YAXPY
#define zzdotu_           ZZDOTU
#define yydotu_           YYDOTU
#define zzdotc_           ZZDOTC
#define yydotc_           YYDOTC
#define dscnrm2_          DSCNRM2
#define dznrm2_           DZNRM2
#define eynrm2_           EYNRM2
#define dscasum_          DSCASUM
#define dzasum_           DZASUM
#define eyasum_           EYASUM
#define izamax_           IZAMAX
#define iyamax_           IYAMAX
                                                     /* Level-2 BLAS */
#define sgemv_            SGEMV
#define egemv_            EGEMV
#define ssymv_            SSYMV
#define esymv_            ESYMV
#define strmv_            STRMV
#define etrmv_            ETRMV
#define strsv_            STRSV
#define etrsv_            ETRSV
#define sger_             SGER
#define eger_             EGER
#define ssyr_             SSYR
#define esyr_             ESYR
#define ssyr2_            SSYR2
#define esyr2_            ESYR2

#define dgemv_            DGEMV
#define egemv_            EGEMV
#define dsymv_            DSYMV
#define esymv_            ESYMV
#define dtrmv_            DTRMV
#define etrmv_            ETRMV
#define dtrsv_            DTRSV
#define etrsv_            ETRSV
#define dger_             DGER
#define eger_             EGER
#define dsyr_             DSYR
#define esyr_             ESYR
#define dsyr2_            DSYR2
#define esyr2_            ESYR2

#define cgemv_            CGEMV
#define ygemv_            YGEMV
#define chemv_            CHEMV
#define yhemv_            YHEMV
#define ctrmv_            CTRMV
#define ytrmv_            YTRMV
#define ctrsv_            CTRSV
#define ytrsv_            YTRSV
#define cgeru_            CGERU
#define ygeru_            YGERU
#define cgerc_            CGERC
#define ygerc_            YGERC
#define cher_             CHER
#define yher_             YHER
#define cher2_            CHER2
#define yher2_            YHER2

#define zgemv_            ZGEMV
#define ygemv_            YGEMV
#define zhemv_            ZHEMV
#define yhemv_            YHEMV
#define ztrmv_            ZTRMV
#define ytrmv_            YTRMV
#define ztrsv_            ZTRSV
#define ytrsv_            YTRSV
#define zgeru_            ZGERU
#define ygeru_            YGERU
#define zgerc_            ZGERC
#define ygerc_            YGERC
#define zher_             ZHER
#define yher_             YHER
#define zher2_            ZHER2
#define yher2_            YHER2
                                                     /* Level-3 BLAS */
#define sgemm_            SGEMM
#define egemm_            EGEMM
#define ssymm_            SSYMM
#define esymm_            ESYMM
#define ssyrk_            SSYRK
#define esyrk_            ESYRK
#define ssyr2k_           SSYR2K
#define esyr2k_           ESYR2K
#define strmm_            STRMM
#define etrmm_            ETRMM
#define strsm_            STRSM
#define etrsm_            ETRSM

#define dgemm_            DGEMM
#define egemm_            EGEMM
#define dsymm_            DSYMM
#define esymm_            ESYMM
#define dsyrk_            DSYRK
#define esyrk_            ESYRK
#define dsyr2k_           DSYR2K
#define esyr2k_           ESYR2K
#define dtrmm_            DTRMM
#define etrmm_            ETRMM
#define dtrsm_            DTRSM
#define etrsm_            ETRSM

#define cgemm_            CGEMM
#define ygemm_            YGEMM
#define chemm_            CHEMM
#define yhemm_            YHEMM
#define csymm_            CSYMM
#define ysymm_            YSYMM
#define csyrk_            CSYRK
#define ysyrk_            YSYRK
#define cherk_            CHERK
#define yherk_            YHERK
#define csyr2k_           CSYR2K
#define ysyr2k_           YSYR2K
#define cher2k_           CHER2K
#define yher2k_           YHER2K
#define ctrmm_            CTRMM
#define ytrmm_            YTRMM
#define ctrsm_            CTRSM
#define ytrsm_            YTRSM

#define zgemm_            ZGEMM
#define ygemm_            YGEMM
#define zhemm_            ZHEMM
#define yhemm_            YHEMM
#define zsymm_            ZSYMM
#define ysymm_            YSYMM
#define zsyrk_            ZSYRK
#define ysyrk_            YSYRK
#define zherk_            ZHERK
#define yherk_            YHERK
#define zsyr2k_           ZSYR2K
#define ysyr2k_           YSYR2K
#define zher2k_           ZHER2K
#define yher2k_           YHER2K
#define ztrmm_            ZTRMM
#define ytrmm_            YTRMM
#define ztrsm_            ZTRSM
#define ytrsm_            YTRSM
                                   /* absolute value auxiliary PBLAS */
#define psatrmv_          PSATRMV
#define peatrmv_          PEATRMV
#define pdatrmv_          PDATRMV
#define peatrmv_          PEATRMV
#define pcatrmv_          PCATRMV
#define pyatrmv_          PYATRMV
#define pzatrmv_          PZATRMV
#define pyatrmv_          PYATRMV
#define psagemv_          PSAGEMV
#define peagemv_          PEAGEMV
#define pdagemv_          PDAGEMV
#define peagemv_          PEAGEMV
#define pcagemv_          PCAGEMV
#define pyagemv_          PYAGEMV
#define pzagemv_          PZAGEMV
#define pyagemv_          PYAGEMV
#define psasymv_          PSASYMV
#define peasymv_          PEASYMV
#define pdasymv_          PDASYMV
#define peasymv_          PEASYMV
#define pcahemv_          PCAHEMV
#define pyahemv_          PYAHEMV
#define pzahemv_          PZAHEMV
#define pyahemv_          PYAHEMV
                                                /* Auxiliary PB-BLAS */
#define pbcmatadd_        PBCMATADD
#define pbymatadd_        PBYMATADD
#define pbdmatadd_        PBDMATADD
#define pbematadd_        PBEMATADD
#define pbsmatadd_        PBSMATADD
#define pbematadd_        PBEMATADD
#define pbzmatadd_        PBZMATADD
#define pbymatadd_        PBYMATADD
                                                   /* Level-2 PBBLAS */
#define pbcgemv_          PBCGEMV
#define pbygemv_          PBYGEMV
#define pbcgeru_          PBCGERU
#define pbygeru_          PBYGERU
#define pbcgerc_          PBCGERC
#define pbygerc_          PBYGERC
#define pbchemv_          PBCHEMV
#define pbyhemv_          PBYHEMV
#define pbcher_           PBCHER
#define pbyher_           PBYHER
#define pbcher2_          PBCHER2
#define pbyher2_          PBYHER2
#define pbctrmv_          PBCTRMV
#define pbytrmv_          PBYTRMV
#define pbctrnv_          PBCTRNV
#define pbytrnv_          PBYTRNV
#define pbctrsv_          PBCTRSV
#define pbytrsv_          PBYTRSV

#define pbdgemv_          PBDGEMV
#define pbegemv_          PBEGEMV
#define pbdger_           PBDGER
#define pbeger_           PBEGER
#define pbdsymv_          PBDSYMV
#define pbesymv_          PBESYMV
#define pbdsyr_           PBDSYR
#define pbesyr_           PBESYR
#define pbdsyr2_          PBDSYR2
#define pbesyr2_          PBESYR2
#define pbdtrmv_          PBDTRMV
#define pbetrmv_          PBETRMV
#define pbdtrnv_          PBDTRNV
#define pbetrnv_          PBETRNV
#define pbdtrsv_          PBDTRSV
#define pbetrsv_          PBETRSV

#define pbsgemv_          PBSGEMV
#define pbegemv_          PBEGEMV
#define pbsger_           PBSGER
#define pbeger_           PBEGER
#define pbssymv_          PBSSYMV
#define pbesymv_          PBESYMV
#define pbssyr_           PBSSYR
#define pbesyr_           PBESYR
#define pbssyr2_          PBSSYR2
#define pbesyr2_          PBESYR2
#define pbstrmv_          PBSTRMV
#define pbetrmv_          PBETRMV
#define pbstrnv_          PBSTRNV
#define pbetrnv_          PBETRNV
#define pbstrsv_          PBSTRSV
#define pbetrsv_          PBETRSV

#define pbzgemv_          PBZGEMV
#define pbygemv_          PBYGEMV
#define pbzgeru_          PBZGERU
#define pbygeru_          PBYGERU
#define pbzgerc_          PBZGERC
#define pbygerc_          PBYGERC
#define pbzhemv_          PBZHEMV
#define pbyhemv_          PBYHEMV
#define pbzher_           PBZHER
#define pbyher_           PBYHER
#define pbzher2_          PBZHER2
#define pbyher2_          PBYHER2
#define pbztrmv_          PBZTRMV
#define pbytrmv_          PBYTRMV
#define pbztrnv_          PBZTRNV
#define pbytrnv_          PBYTRNV
#define pbztrsv_          PBZTRSV
#define pbytrsv_          PBYTRSV
                                                   /* Level-3 PBBLAS */
#define pbcgemm_          PBCGEMM
#define pbygemm_          PBYGEMM
#define pbchemm_          PBCHEMM
#define pbyhemm_          PBYHEMM
#define pbcher2k_         PBCHER2K
#define pbyher2k_         PBYHER2K
#define pbcherk_          PBCHERK
#define pbyherk_          PBYHERK
#define pbcsymm_          PBCSYMM
#define pbysymm_          PBYSYMM
#define pbcsyr2k_         PBCSYR2K
#define pbysyr2k_         PBYSYR2K
#define pbcsyrk_          PBCSYRK
#define pbysyrk_          PBYSYRK
#define pbctrmm_          PBCTRMM
#define pbytrmm_          PBYTRMM
#define pbctrsm_          PBCTRSM
#define pbytrsm_          PBYTRSM
#define pbctran_          PBCTRAN
#define pbytran_          PBYTRAN

#define pbdgemm_          PBDGEMM
#define pbegemm_          PBEGEMM
#define pbdsymm_          PBDSYMM
#define pbesymm_          PBESYMM
#define pbdsyr2k_         PBDSYR2K
#define pbesyr2k_         PBESYR2K
#define pbdsyrk_          PBDSYRK
#define pbesyrk_          PBESYRK
#define pbdtrmm_          PBDTRMM
#define pbetrmm_          PBETRMM
#define pbdtrsm_          PBDTRSM
#define pbetrsm_          PBETRSM
#define pbdtran_          PBDTRAN
#define pbetran_          PBETRAN

#define pbsgemm_          PBSGEMM
#define pbegemm_          PBEGEMM
#define pbssymm_          PBSSYMM
#define pbesymm_          PBESYMM
#define pbssyr2k_         PBSSYR2K
#define pbesyr2k_         PBESYR2K
#define pbssyrk_          PBSSYRK
#define pbesyrk_          PBESYRK
#define pbstrmm_          PBSTRMM
#define pbetrmm_          PBETRMM
#define pbstrsm_          PBSTRSM
#define pbetrsm_          PBETRSM
#define pbstran_          PBSTRAN
#define pbetran_          PBETRAN

#define pbzgemm_          PBZGEMM
#define pbygemm_          PBYGEMM
#define pbzhemm_          PBZHEMM
#define pbyhemm_          PBYHEMM
#define pbzher2k_         PBZHER2K
#define pbyher2k_         PBYHER2K
#define pbzherk_          PBZHERK
#define pbyherk_          PBYHERK
#define pbzsymm_          PBZSYMM
#define pbysymm_          PBYSYMM
#define pbzsyr2k_         PBZSYR2K
#define pbysyr2k_         PBYSYR2K
#define pbzsyrk_          PBZSYRK
#define pbysyrk_          PBYSYRK
#define pbztrmm_          PBZTRMM
#define pbytrmm_          PBYTRMM
#define pbztrsm_          PBZTRSM
#define pbytrsm_          PBYTRSM
#define pbztran_          PBZTRAN
#define pbytran_          PBYTRAN
                                                 /* Auxilliary PBLAS */
#define pberror_          PBERROR
#define pb_freebuf_       PB_FREEBUF
#define pb_topget_        PB_TOPGET
#define pb_topset_        PB_TOPSET
                                                    /* Level-1 PBLAS */
#define psrotg_           PSROTG
#define perotg_           PEROTG
#define psrotmg_          PSROTMG
#define perotmg_          PEROTMG
#define psrot_            PSROT
#define perot_            PEROT
#define psrotm_           PSROTM
#define perotm_           PEROTM
#define psswap_           PSSWAP
#define peswap_           PESWAP
#define psscal_           PSSCAL
#define pescal_           PESCAL
#define pscopy_           PSCOPY
#define pecopy_           PECOPY
#define psaxpy_           PSAXPY
#define peaxpy_           PEAXPY
#define psdot_            PSDOT
#define pedot_            PEDOT
#define psnrm2_           PSNRM2
#define penrm2_           PENRM2
#define psasum_           PSASUM
#define peasum_           PEASUM
#define psamax_           PSAMAX
#define peamax_           PEAMAX

#define pdrotg_           PDROTG
#define perotg_           PEROTG
#define pdrotmg_          PDROTMG
#define perotmg_          PEROTMG
#define pdrot_            PDROT
#define perot_            PEROT
#define pdrotm_           PDROTM
#define perotm_           PEROTM
#define pdswap_           PDSWAP
#define peswap_           PESWAP
#define pdscal_           PDSCAL
#define pescal_           PESCAL
#define pdcopy_           PDCOPY
#define pecopy_           PECOPY
#define pdaxpy_           PDAXPY
#define peaxpy_           PEAXPY
#define pddot_            PDDOT
#define pedot_            PEDOT
#define pdnrm2_           PDNRM2
#define penrm2_           PENRM2
#define pdasum_           PDASUM
#define peasum_           PEASUM
#define pdamax_           PDAMAX
#define peamax_           PEAMAX

#define pcswap_           PCSWAP
#define pyswap_           PYSWAP
#define pcscal_           PCSCAL
#define pyscal_           PYSCAL
#define pcsscal_          PCSSCAL
#define pyescal_          PYESCAL
#define pccopy_           PCCOPY
#define pycopy_           PYCOPY
#define pcaxpy_           PCAXPY
#define pyaxpy_           PYAXPY
#define pcdotu_           PCDOTU
#define pydotu_           PYDOTU
#define pcdotc_           PCDOTC
#define pydotc_           PYDOTC
#define pscnrm2_          PSCNRM2
#define peynrm2_          PEYNRM2
#define pscasum_          PSCASUM
#define peyasum_          PEYASUM
#define pcamax_           PCAMAX
#define pyamax_           PYAMAX
#define pcrot_            PCROT
#define pyrot_            PYROT
#define crot_             CROT
#define yrot_             YROT

#define pzswap_           PZSWAP
#define pyswap_           PYSWAP
#define pzscal_           PZSCAL
#define pyscal_           PYSCAL
#define pzdscal_          PZDSCAL
#define pyescal_          PYESCAL
#define pzcopy_           PZCOPY
#define pycopy_           PYCOPY
#define pzaxpy_           PZAXPY
#define pyaxpy_           PYAXPY
#define pzdotu_           PZDOTU
#define pydotu_           PYDOTU
#define pzdotc_           PZDOTC
#define pydotc_           PYDOTC
#define pdznrm2_          PDZNRM2
#define peynrm2_          PEYNRM2
#define pdzasum_          PDZASUM
#define peyasum_          PEYASUM
#define pzamax_           PZAMAX
#define pyamax_           PYAMAX
#define pzrot_            PZROT
#define pyrot_            PYROT
#define zrot_             ZROT
#define yrot_             YROT
                                                    /* Level-2 PBLAS */
#define pcgemv_           PCGEMV
#define pygemv_           PYGEMV
#define pcgeru_           PCGERU
#define pygeru_           PYGERU
#define pcgerc_           PCGERC
#define pygerc_           PYGERC
#define pchemv_           PCHEMV
#define pyhemv_           PYHEMV
#define pcher_            PCHER
#define pyher_            PYHER
#define pcher2_           PCHER2
#define pyher2_           PYHER2
#define pctrmv_           PCTRMV
#define pytrmv_           PYTRMV
#define pctrsv_           PCTRSV
#define pytrsv_           PYTRSV

#define pdgemv_           PDGEMV
#define pegemv_           PEGEMV
#define pdger_            PDGER
#define peger_            PEGER
#define pdsymv_           PDSYMV
#define pesymv_           PESYMV
#define pdsyr_            PDSYR
#define pesyr_            PESYR
#define pdsyr2_           PDSYR2
#define pesyr2_           PESYR2
#define pdtrmv_           PDTRMV
#define petrmv_           PETRMV
#define pdtrsv_           PDTRSV
#define petrsv_           PETRSV

#define psgemv_           PSGEMV
#define pegemv_           PEGEMV
#define psger_            PSGER
#define peger_            PEGER
#define pssymv_           PSSYMV
#define pesymv_           PESYMV
#define pssyr_            PSSYR
#define pesyr_            PESYR
#define pssyr2_           PSSYR2
#define pesyr2_           PESYR2
#define pstrmv_           PSTRMV
#define petrmv_           PETRMV
#define pstrsv_           PSTRSV
#define petrsv_           PETRSV

#define pzgemv_           PZGEMV
#define pygemv_           PYGEMV
#define pzgeru_           PZGERU
#define pygeru_           PYGERU
#define pzgerc_           PZGERC
#define pygerc_           PYGERC
#define pzhemv_           PZHEMV
#define pyhemv_           PYHEMV
#define pzher_            PZHER
#define pyher_            PYHER
#define pzher2_           PZHER2
#define pyher2_           PYHER2
#define pztrmv_           PZTRMV
#define pytrmv_           PYTRMV
#define pztrsv_           PZTRSV
#define pytrsv_           PYTRSV
                                                    /* Level-3 PBLAS */
#define pcgemm_           PCGEMM
#define pygemm_           PYGEMM
#define pchemm_           PCHEMM
#define pyhemm_           PYHEMM
#define pcher2k_          PCHER2K
#define pyher2k_          PYHER2K
#define pcherk_           PCHERK
#define pyherk_           PYHERK
#define pcsymm_           PCSYMM
#define pysymm_           PYSYMM
#define pcsyr2k_          PCSYR2K
#define pysyr2k_          PYSYR2K
#define pcsyrk_           PCSYRK
#define pysyrk_           PYSYRK
#define pctrmm_           PCTRMM
#define pytrmm_           PYTRMM
#define pctrsm_           PCTRSM
#define pytrsm_           PYTRSM
#define pctranu_          PCTRANU
#define pytranu_          PYTRANU
#define pctranc_          PCTRANC
#define pytranc_          PYTRANC

#define pdgemm_           PDGEMM
#define pegemm_           PEGEMM
#define pdsymm_           PDSYMM
#define pesymm_           PESYMM
#define pdsyr2k_          PDSYR2K
#define pesyr2k_          PESYR2K
#define pdsyrk_           PDSYRK
#define pesyrk_           PESYRK
#define pdtrmm_           PDTRMM
#define petrmm_           PETRMM
#define pdtrsm_           PDTRSM
#define petrsm_           PETRSM
#define pdtran_           PDTRAN
#define petran_           PETRAN

#define psgemm_           PSGEMM
#define pegemm_           PEGEMM
#define pssymm_           PSSYMM
#define pesymm_           PESYMM
#define pssyr2k_          PSSYR2K
#define pesyr2k_          PESYR2K
#define pssyrk_           PSSYRK
#define pesyrk_           PESYRK
#define pstrmm_           PSTRMM
#define petrmm_           PETRMM
#define pstrsm_           PSTRSM
#define petrsm_           PETRSM
#define pstran_           PSTRAN
#define petran_           PETRAN

#define pzgemm_           PZGEMM
#define pygemm_           PYGEMM
#define pzhemm_           PZHEMM
#define pyhemm_           PYHEMM
#define pzher2k_          PZHER2K
#define pyher2k_          PYHER2K
#define pzherk_           PZHERK
#define pyherk_           PYHERK
#define pzsymm_           PZSYMM
#define pysymm_           PYSYMM
#define pzsyr2k_          PZSYR2K
#define pysyr2k_          PYSYR2K
#define pzsyrk_           PZSYRK
#define pysyrk_           PYSYRK
#define pztrmm_           PZTRMM
#define pytrmm_           PYTRMM
#define pztrsm_           PZTRSM
#define pytrsm_           PYTRSM
#define pztranu_          PZTRANU
#define pytranu_          PYTRANU
#define pztranc_          PZTRANC
#define pytranc_          PYTRANC

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* for following FORTRAN to C interface:
*           FORTRAN CALL               C DECLARATION
*           call pdgemm(...)           void pdgemm(...)
*           call pegemm(...)           void pegemm(...)
*/
                                                            /* TOOLS */
#define ilcm_             ilcm
#define infog2l_          infog2l
#define numroc_           numroc
#define pstreecomb_       pstreecomb
#define petreecomb_       petreecomb
#define pdtreecomb_       pdtreecomb
#define petreecomb_       petreecomb
#define pctreecomb_       pctreecomb
#define pytreecomb_       pytreecomb
#define pztreecomb_       pztreecomb
#define pytreecomb_       pytreecomb
#define scombamax_        scombamax
#define ecombamax_        ecombamax
#define dcombamax_        dcombamax
#define ecombamax_        ecombamax
#define ccombamax_        ccombamax
#define ycombamax_        ycombamax
#define zcombamax_        zcombamax
#define ycombamax_        ycombamax
#define scombnrm2_        scombnrm2
#define ecombnrm2_        ecombnrm2
#define dcombnrm2_        dcombnrm2
#define ecombnrm2_        ecombnrm2

#define dlamov_           dlamov
#define elamov_           elamov
#define slamov_           slamov
#define elamov_           elamov
#define clamov_           clamov
#define ylamov_           ylamov
#define zlamov_           zlamov
#define ylamov_           ylamov
#define dlacpy_           dlacpy
#define elacpy_           elacpy
#define slacpy_           slacpy
#define elacpy_           elacpy
#define clacpy_           clacpy
#define ylacpy_           ylacpy
#define zlacpy_           zlacpy
#define ylacpy_           ylacpy
#define xerbla_           xerbla
                                                            /* BLACS */
#define blacs_abort_      blacs_abort
#define blacs_gridinfo_   blacs_gridinfo

#define igesd2d_          igesd2d
#define igebs2d_          igebs2d
#define itrsd2d_          itrsd2d
#define itrbs2d_          itrbs2d
#define igerv2d_          igerv2d
#define igebr2d_          igebr2d
#define itrrv2d_          itrrv2d
#define itrbr2d_          itrbr2d
#define igamx2d_          igamx2d
#define igamn2d_          igamn2d
#define igsum2d_          igsum2d

#define sgesd2d_          sgesd2d
#define egesd2d_          egesd2d
#define sgebs2d_          sgebs2d
#define egebs2d_          egebs2d
#define strsd2d_          strsd2d
#define etrsd2d_          etrsd2d
#define strbs2d_          strbs2d
#define etrbs2d_          etrbs2d
#define sgerv2d_          sgerv2d
#define egerv2d_          egerv2d
#define sgebr2d_          sgebr2d
#define egebr2d_          egebr2d
#define strrv2d_          strrv2d
#define etrrv2d_          etrrv2d
#define strbr2d_          strbr2d
#define etrbr2d_          etrbr2d
#define sgamx2d_          sgamx2d
#define egamx2d_          egamx2d
#define sgamn2d_          sgamn2d
#define egamn2d_          egamn2d
#define sgsum2d_          sgsum2d
#define egsum2d_          egsum2d

#define dgesd2d_          dgesd2d
#define egesd2d_          egesd2d
#define dgebs2d_          dgebs2d
#define egebs2d_          egebs2d
#define dtrsd2d_          dtrsd2d
#define etrsd2d_          etrsd2d
#define dtrbs2d_          dtrbs2d
#define etrbs2d_          etrbs2d
#define dgerv2d_          dgerv2d
#define egerv2d_          egerv2d
#define dgebr2d_          dgebr2d
#define egebr2d_          egebr2d
#define dtrrv2d_          dtrrv2d
#define etrrv2d_          etrrv2d
#define dtrbr2d_          dtrbr2d
#define etrbr2d_          etrbr2d
#define dgamx2d_          dgamx2d
#define egamx2d_          egamx2d
#define dgamn2d_          dgamn2d
#define egamn2d_          egamn2d
#define dgsum2d_          dgsum2d
#define egsum2d_          egsum2d

#define cgesd2d_          cgesd2d
#define ygesd2d_          ygesd2d
#define cgebs2d_          cgebs2d
#define ygebs2d_          ygebs2d
#define ctrsd2d_          ctrsd2d
#define ytrsd2d_          ytrsd2d
#define ctrbs2d_          ctrbs2d
#define ytrbs2d_          ytrbs2d
#define cgerv2d_          cgerv2d
#define ygerv2d_          ygerv2d
#define cgebr2d_          cgebr2d
#define ygebr2d_          ygebr2d
#define ctrrv2d_          ctrrv2d
#define ytrrv2d_          ytrrv2d
#define ctrbr2d_          ctrbr2d
#define ytrbr2d_          ytrbr2d
#define cgamx2d_          cgamx2d
#define ygamx2d_          ygamx2d
#define cgamn2d_          cgamn2d
#define ygamn2d_          ygamn2d
#define cgsum2d_          cgsum2d
#define ygsum2d_          ygsum2d

#define zgesd2d_          zgesd2d
#define ygesd2d_          ygesd2d
#define zgebs2d_          zgebs2d
#define ygebs2d_          ygebs2d
#define ztrsd2d_          ztrsd2d
#define ytrsd2d_          ytrsd2d
#define ztrbs2d_          ztrbs2d
#define ytrbs2d_          ytrbs2d
#define zgerv2d_          zgerv2d
#define ygerv2d_          ygerv2d
#define zgebr2d_          zgebr2d
#define ygebr2d_          ygebr2d
#define ztrrv2d_          ztrrv2d
#define ytrrv2d_          ytrrv2d
#define ztrbr2d_          ztrbr2d
#define ytrbr2d_          ytrbr2d
#define zgamx2d_          zgamx2d
#define ygamx2d_          ygamx2d
#define zgamn2d_          zgamn2d
#define ygamn2d_          ygamn2d
#define zgsum2d_          zgsum2d
#define ygsum2d_          ygsum2d
                                                     /* Level-1 BLAS */
#define srotg_            srotg
#define erotg_            erotg
#define srotmg_           srotmg
#define erotmg_           erotmg
#define srot_             srot
#define erot_             erot
#define srotm_            srotm
#define erotm_            erotm
#define sswap_            sswap
#define eswap_            eswap
#define sscal_            sscal
#define escal_            escal
#define scopy_            scopy
#define ecopy_            ecopy
#define saxpy_            saxpy
#define eaxpy_            eaxpy
#define ssdot_            ssdot
#define esdot_            esdot
#define isamax_           isamax
#define ieamax_           ieamax

#define drotg_            drotg
#define erotg_            erotg
#define drotmg_           drotmg
#define erotmg_           erotmg
#define drot_             drot
#define erot_             erot
#define drotm_            drotm
#define erotm_            erotm
#define dswap_            dswap
#define eswap_            eswap
#define dscal_            dscal
#define escal_            escal
#define dcopy_            dcopy
#define ecopy_            ecopy
#define daxpy_            daxpy
#define eaxpy_            eaxpy
#define dddot_            dddot
#define dnrm2_            dnrm2
#define enrm2_            enrm2
#define dsnrm2_           dsnrm2
#define dasum_            dasum
#define easum_            easum
#define dsasum_           dsasum
#define idamax_           idamax
#define ieamax_           ieamax

#define cswap_            cswap
#define yswap_            yswap
#define cscal_            cscal
#define yscal_            yscal
#define csscal_           csscal
#define yescal_           yescal
#define ccopy_            ccopy
#define ycopy_            ycopy
#define caxpy_            caxpy
#define yaxpy_            yaxpy
#define ccdotu_           ccdotu
#define yydotu_           yydotu
#define ccdotc_           ccdotc
#define yydotc_           yydotc
#define icamax_           icamax
#define iyamax_           iyamax

#define zswap_            zswap
#define yswap_            yswap
#define zscal_            zscal
#define yscal_            yscal
#define zdscal_           zdscal
#define yescal_           yescal
#define zcopy_            zcopy
#define ycopy_            ycopy
#define zaxpy_            zaxpy
#define yaxpy_            yaxpy
#define zzdotu_           zzdotu
#define yydotu_           yydotu
#define zzdotc_           zzdotc
#define yydotc_           yydotc
#define dscnrm2_          dscnrm2
#define dznrm2_           dznrm2
#define eynrm2_           eynrm2
#define dscasum_          dscasum
#define dzasum_           dzasum
#define eyasum_           eyasum
#define izamax_           izamax
#define iyamax_           iyamax
                                                     /* Level-2 BLAS */
#define sgemv_            sgemv
#define egemv_            egemv
#define ssymv_            ssymv
#define esymv_            esymv
#define strmv_            strmv
#define etrmv_            etrmv
#define strsv_            strsv
#define etrsv_            etrsv
#define sger_             sger
#define eger_             eger
#define ssyr_             ssyr
#define esyr_             esyr
#define ssyr2_            ssyr2
#define esyr2_            esyr2

#define dgemv_            dgemv
#define egemv_            egemv
#define dsymv_            dsymv
#define esymv_            esymv
#define dtrmv_            dtrmv
#define etrmv_            etrmv
#define dtrsv_            dtrsv
#define etrsv_            etrsv
#define dger_             dger
#define eger_             eger
#define dsyr_             dsyr
#define esyr_             esyr
#define dsyr2_            dsyr2
#define esyr2_            esyr2

#define cgemv_            cgemv
#define ygemv_            ygemv
#define chemv_            chemv
#define yhemv_            yhemv
#define ctrmv_            ctrmv
#define ytrmv_            ytrmv
#define ctrsv_            ctrsv
#define ytrsv_            ytrsv
#define cgeru_            cgeru
#define ygeru_            ygeru
#define cgerc_            cgerc
#define ygerc_            ygerc
#define cher_             cher
#define yher_             yher
#define cher2_            cher2
#define yher2_            yher2

#define zgemv_            zgemv
#define ygemv_            ygemv
#define zhemv_            zhemv
#define yhemv_            yhemv
#define ztrmv_            ztrmv
#define ytrmv_            ytrmv
#define ztrsv_            ztrsv
#define ytrsv_            ytrsv
#define zgeru_            zgeru
#define ygeru_            ygeru
#define zgerc_            zgerc
#define ygerc_            ygerc
#define zher_             zher
#define yher_             yher
#define zher2_            zher2
#define yher2_            yher2
                                                     /* Level-3 BLAS */
#define sgemm_            sgemm
#define egemm_            egemm
#define ssymm_            ssymm
#define esymm_            esymm
#define ssyrk_            ssyrk
#define esyrk_            esyrk
#define ssyr2k_           ssyr2k
#define esyr2k_           esyr2k
#define strmm_            strmm
#define etrmm_            etrmm
#define strsm_            strsm
#define etrsm_            etrsm

#define dgemm_            dgemm
#define egemm_            egemm
#define dsymm_            dsymm
#define esymm_            esymm
#define dsyrk_            dsyrk
#define esyrk_            esyrk
#define dsyr2k_           dsyr2k
#define esyr2k_           esyr2k
#define dtrmm_            dtrmm
#define etrmm_            etrmm
#define dtrsm_            dtrsm
#define etrsm_            etrsm

#define cgemm_            cgemm
#define ygemm_            ygemm
#define chemm_            chemm
#define yhemm_            yhemm
#define csymm_            csymm
#define ysymm_            ysymm
#define csyrk_            csyrk
#define ysyrk_            ysyrk
#define cherk_            cherk
#define yherk_            yherk
#define csyr2k_           csyr2k
#define ysyr2k_           ysyr2k
#define cher2k_           cher2k
#define yher2k_           yher2k
#define ctrmm_            ctrmm
#define ytrmm_            ytrmm
#define ctrsm_            ctrsm
#define ytrsm_            ytrsm

#define zgemm_            zgemm
#define ygemm_            ygemm
#define zhemm_            zhemm
#define yhemm_            yhemm
#define zsymm_            zsymm
#define ysymm_            ysymm
#define zsyrk_            zsyrk
#define ysyrk_            ysyrk
#define zherk_            zherk
#define yherk_            yherk
#define zsyr2k_           zsyr2k
#define ysyr2k_           ysyr2k
#define zher2k_           zher2k
#define yher2k_           yher2k
#define ztrmm_            ztrmm
#define ytrmm_            ytrmm
#define ztrsm_            ztrsm
#define ytrsm_            ytrsm
                                   /* absolute value auxiliary PBLAS */
#define psatrmv_          psatrmv
#define peatrmv_          peatrmv
#define pdatrmv_          pdatrmv
#define peatrmv_          peatrmv
#define pcatrmv_          pcatrmv
#define pyatrmv_          pyatrmv
#define pzatrmv_          pzatrmv
#define pyatrmv_          pyatrmv
#define psagemv_          psagemv
#define peagemv_          peagemv
#define pdagemv_          pdagemv
#define peagemv_          peagemv
#define pcagemv_          pcagemv
#define pyagemv_          pyagemv
#define pzagemv_          pzagemv
#define pyagemv_          pyagemv
#define psasymv_          psasymv
#define peasymv_          peasymv
#define pdasymv_          pdasymv
#define peasymv_          peasymv
#define pcahemv_          pcahemv
#define pyahemv_          pyahemv
#define pzahemv_          pzahemv
#define pyahemv_          pyahemv
                                                /* Auxiliary PB-BLAS */
#define pbcmatadd_        pbcmatadd
#define pbymatadd_        pbymatadd
#define pbdmatadd_        pbdmatadd
#define pbematadd_        pbematadd
#define pbsmatadd_        pbsmatadd
#define pbematadd_        pbematadd
#define pbzmatadd_        pbzmatadd
#define pbymatadd_        pbymatadd
                                                   /* Level-2 PBBLAS */
#define pbcgemv_          pbcgemv
#define pbygemv_          pbygemv
#define pbcgeru_          pbcgeru
#define pbygeru_          pbygeru
#define pbcgerc_          pbcgerc
#define pbygerc_          pbygerc
#define pbchemv_          pbchemv
#define pbyhemv_          pbyhemv
#define pbcher_           pbcher
#define pbyher_           pbyher
#define pbcher2_          pbcher2
#define pbyher2_          pbyher2
#define pbctrmv_          pbctrmv
#define pbytrmv_          pbytrmv
#define pbctrnv_          pbctrnv
#define pbytrnv_          pbytrnv
#define pbctrsv_          pbctrsv
#define pbytrsv_          pbytrsv

#define pbdgemv_          pbdgemv
#define pbegemv_          pbegemv
#define pbdger_           pbdger
#define pbeger_           pbeger
#define pbdsymv_          pbdsymv
#define pbesymv_          pbesymv
#define pbdsyr_           pbdsyr
#define pbesyr_           pbesyr
#define pbdsyr2_          pbdsyr2
#define pbesyr2_          pbesyr2
#define pbdtrmv_          pbdtrmv
#define pbetrmv_          pbetrmv
#define pbdtrnv_          pbdtrnv
#define pbetrnv_          pbetrnv
#define pbdtrsv_          pbdtrsv
#define pbetrsv_          pbetrsv

#define pbsgemv_          pbsgemv
#define pbegemv_          pbegemv
#define pbsger_           pbsger
#define pbeger_           pbeger
#define pbssymv_          pbssymv
#define pbesymv_          pbesymv
#define pbssyr_           pbssyr
#define pbesyr_           pbesyr
#define pbssyr2_          pbssyr2
#define pbesyr2_          pbesyr2
#define pbstrmv_          pbstrmv
#define pbetrmv_          pbetrmv
#define pbstrnv_          pbstrnv
#define pbetrnv_          pbetrnv
#define pbstrsv_          pbstrsv
#define pbetrsv_          pbetrsv

#define pbzgemv_          pbzgemv
#define pbygemv_          pbygemv
#define pbzgeru_          pbzgeru
#define pbygeru_          pbygeru
#define pbzgerc_          pbzgerc
#define pbygerc_          pbygerc
#define pbzhemv_          pbzhemv
#define pbyhemv_          pbyhemv
#define pbzher_           pbzher
#define pbyher_           pbyher
#define pbzher2_          pbzher2
#define pbyher2_          pbyher2
#define pbztrmv_          pbztrmv
#define pbytrmv_          pbytrmv
#define pbztrnv_          pbztrnv
#define pbytrnv_          pbytrnv
#define pbztrsv_          pbztrsv
#define pbytrsv_          pbytrsv
                                                   /* Level-3 PBBLAS */
#define pbcgemm_          pbcgemm
#define pbygemm_          pbygemm
#define pbchemm_          pbchemm
#define pbyhemm_          pbyhemm
#define pbcher2k_         pbcher2k
#define pbyher2k_         pbyher2k
#define pbcherk_          pbcherk
#define pbyherk_          pbyherk
#define pbcsymm_          pbcsymm
#define pbysymm_          pbysymm
#define pbcsyr2k_         pbcsyr2k
#define pbysyr2k_         pbysyr2k
#define pbcsyrk_          pbcsyrk
#define pbysyrk_          pbysyrk
#define pbctrmm_          pbctrmm
#define pbytrmm_          pbytrmm
#define pbctrsm_          pbctrsm
#define pbytrsm_          pbytrsm
#define pbctran_          pbctran
#define pbytran_          pbytran

#define pbdgemm_          pbdgemm
#define pbegemm_          pbegemm
#define pbdsymm_          pbdsymm
#define pbesymm_          pbesymm
#define pbdsyr2k_         pbdsyr2k
#define pbesyr2k_         pbesyr2k
#define pbdsyrk_          pbdsyrk
#define pbesyrk_          pbesyrk
#define pbdtrmm_          pbdtrmm
#define pbetrmm_          pbetrmm
#define pbdtrsm_          pbdtrsm
#define pbetrsm_          pbetrsm
#define pbdtran_          pbdtran
#define pbetran_          pbetran

#define pbsgemm_          pbsgemm
#define pbegemm_          pbegemm
#define pbssymm_          pbssymm
#define pbesymm_          pbesymm
#define pbssyr2k_         pbssyr2k
#define pbesyr2k_         pbesyr2k
#define pbssyrk_          pbssyrk
#define pbesyrk_          pbesyrk
#define pbstrmm_          pbstrmm
#define pbetrmm_          pbetrmm
#define pbstrsm_          pbstrsm
#define pbetrsm_          pbetrsm
#define pbstran_          pbstran
#define pbetran_          pbetran

#define pbzgemm_          pbzgemm
#define pbygemm_          pbygemm
#define pbzhemm_          pbzhemm
#define pbyhemm_          pbyhemm
#define pbzher2k_         pbzher2k
#define pbyher2k_         pbyher2k
#define pbzherk_          pbzherk
#define pbyherk_          pbyherk
#define pbzsymm_          pbzsymm
#define pbysymm_          pbysymm
#define pbzsyr2k_         pbzsyr2k
#define pbysyr2k_         pbysyr2k
#define pbzsyrk_          pbzsyrk
#define pbysyrk_          pbysyrk
#define pbztrmm_          pbztrmm
#define pbytrmm_          pbytrmm
#define pbztrsm_          pbztrsm
#define pbytrsm_          pbytrsm
#define pbztran_          pbztran
#define pbytran_          pbytran
                                                 /* Auxilliary PBLAS */
#define pberror_          pberror
#define pb_freebuf_       pb_freebuf
#define pb_topget_        pb_topget
#define pb_topset_        pb_topset
                                                    /* Level-1 PBLAS */
#define psrotg_           psrotg
#define perotg_           perotg
#define psrotmg_          psrotmg
#define perotmg_          perotmg
#define psrot_            psrot
#define perot_            perot
#define psrotm_           psrotm
#define perotm_           perotm
#define psswap_           psswap
#define peswap_           peswap
#define psscal_           psscal
#define pescal_           pescal
#define pscopy_           pscopy
#define pecopy_           pecopy
#define psaxpy_           psaxpy
#define peaxpy_           peaxpy
#define psdot_            psdot
#define pedot_            pedot
#define psnrm2_           psnrm2
#define penrm2_           penrm2
#define psasum_           psasum
#define peasum_           peasum
#define psamax_           psamax
#define peamax_           peamax

#define pdrotg_           pdrotg
#define perotg_           perotg
#define pdrotmg_          pdrotmg
#define perotmg_          perotmg
#define pdrot_            pdrot
#define perot_            perot
#define pdrotm_           pdrotm
#define perotm_           perotm
#define pdswap_           pdswap
#define peswap_           peswap
#define pdscal_           pdscal
#define pescal_           pescal
#define pdcopy_           pdcopy
#define pecopy_           pecopy
#define pdaxpy_           pdaxpy
#define peaxpy_           peaxpy
#define pddot_            pddot
#define pedot_            pedot
#define pdnrm2_           pdnrm2
#define penrm2_           penrm2
#define pdasum_           pdasum
#define peasum_           peasum
#define pdamax_           pdamax
#define peamax_           peamax

#define pcswap_           pcswap
#define pyswap_           pyswap
#define pcscal_           pcscal
#define pyscal_           pyscal
#define pcsscal_          pcsscal
#define pyescal_          pyescal
#define pccopy_           pccopy
#define pycopy_           pycopy
#define pcaxpy_           pcaxpy
#define pyaxpy_           pyaxpy
#define pcdotu_           pcdotu
#define pydotu_           pydotu
#define pcdotc_           pcdotc
#define pydotc_           pydotc
#define pscnrm2_          pscnrm2
#define peynrm2_          peynrm2
#define pscasum_          pscasum
#define peyasum_          peyasum
#define pcamax_           pcamax
#define pyamax_           pyamax
#define pcrot_            pcrot
#define pyrot_            pyrot
#define crot_             crot
#define yrot_             yrot

#define pzswap_           pzswap
#define pyswap_           pyswap
#define pzscal_           pzscal
#define pyscal_           pyscal
#define pzdscal_          pzdscal
#define pyescal_          pyescal
#define pzcopy_           pzcopy
#define pycopy_           pycopy
#define pzaxpy_           pzaxpy
#define pyaxpy_           pyaxpy
#define pzdotu_           pzdotu
#define pydotu_           pydotu
#define pzdotc_           pzdotc
#define pydotc_           pydotc
#define pdznrm2_          pdznrm2
#define peynrm2_          peynrm2
#define pdzasum_          pdzasum
#define peyasum_          peyasum
#define pzamax_           pzamax
#define pyamax_           pyamax
#define pzrot_            pzrot
#define pyrot_            pyrot
#define zrot_             zrot
#define yrot_             yrot
                                                    /* Level-2 PBLAS */
#define pcgemv_           pcgemv
#define pygemv_           pygemv
#define pcgeru_           pcgeru
#define pygeru_           pygeru
#define pcgerc_           pcgerc
#define pygerc_           pygerc
#define pchemv_           pchemv
#define pyhemv_           pyhemv
#define pcher_            pcher
#define pyher_            pyher
#define pcher2_           pcher2
#define pyher2_           pyher2
#define pctrmv_           pctrmv
#define pytrmv_           pytrmv
#define pctrsv_           pctrsv
#define pytrsv_           pytrsv

#define pdgemv_           pdgemv
#define pegemv_           pegemv
#define pdger_            pdger
#define peger_            peger
#define pdsymv_           pdsymv
#define pesymv_           pesymv
#define pdsyr_            pdsyr
#define pesyr_            pesyr
#define pdsyr2_           pdsyr2
#define pesyr2_           pesyr2
#define pdtrmv_           pdtrmv
#define petrmv_           petrmv
#define pdtrsv_           pdtrsv
#define petrsv_           petrsv

#define psgemv_           psgemv
#define pegemv_           pegemv
#define psger_            psger
#define peger_            peger
#define pssymv_           pssymv
#define pesymv_           pesymv
#define pssyr_            pssyr
#define pesyr_            pesyr
#define pssyr2_           pssyr2
#define pesyr2_           pesyr2
#define pstrmv_           pstrmv
#define petrmv_           petrmv
#define pstrsv_           pstrsv
#define petrsv_           petrsv

#define pzgemv_           pzgemv
#define pygemv_           pygemv
#define pzgeru_           pzgeru
#define pygeru_           pygeru
#define pzgerc_           pzgerc
#define pygerc_           pygerc
#define pzhemv_           pzhemv
#define pyhemv_           pyhemv
#define pzher_            pzher
#define pyher_            pyher
#define pzher2_           pzher2
#define pyher2_           pyher2
#define pztrmv_           pztrmv
#define pytrmv_           pytrmv
#define pztrsv_           pztrsv
#define pytrsv_           pytrsv
                                                    /* Level-3 PBLAS */
#define pcgemm_           pcgemm
#define pygemm_           pygemm
#define pchemm_           pchemm
#define pyhemm_           pyhemm
#define pcher2k_          pcher2k
#define pyher2k_          pyher2k
#define pcherk_           pcherk
#define pyherk_           pyherk
#define pcsymm_           pcsymm
#define pysymm_           pysymm
#define pcsyr2k_          pcsyr2k
#define pysyr2k_          pysyr2k
#define pcsyrk_           pcsyrk
#define pysyrk_           pysyrk
#define pctrmm_           pctrmm
#define pytrmm_           pytrmm
#define pctrsm_           pctrsm
#define pytrsm_           pytrsm
#define pctranu_          pctranu
#define pytranu_          pytranu
#define pctranc_          pctranc
#define pytranc_          pytranc

#define pdgemm_           pdgemm
#define pegemm_           pegemm
#define pdsymm_           pdsymm
#define pesymm_           pesymm
#define pdsyr2k_          pdsyr2k
#define pesyr2k_          pesyr2k
#define pdsyrk_           pdsyrk
#define pesyrk_           pesyrk
#define pdtrmm_           pdtrmm
#define petrmm_           petrmm
#define pdtrsm_           pdtrsm
#define petrsm_           petrsm
#define pdtran_           pdtran
#define petran_           petran

#define psgemm_           psgemm
#define pegemm_           pegemm
#define pssymm_           pssymm
#define pesymm_           pesymm
#define pssyr2k_          pssyr2k
#define pesyr2k_          pesyr2k
#define pssyrk_           pssyrk
#define pesyrk_           pesyrk
#define pstrmm_           pstrmm
#define petrmm_           petrmm
#define pstrsm_           pstrsm
#define petrsm_           petrsm
#define pstran_           pstran
#define petran_           petran

#define pzgemm_           pzgemm
#define pygemm_           pygemm
#define pzhemm_           pzhemm
#define pyhemm_           pyhemm
#define pzher2k_          pzher2k
#define pyher2k_          pyher2k
#define pzherk_           pzherk
#define pyherk_           pyherk
#define pzsymm_           pzsymm
#define pysymm_           pysymm
#define pzsyr2k_          pzsyr2k
#define pysyr2k_          pysyr2k
#define pzsyrk_           pzsyrk
#define pysyrk_           pysyrk
#define pztrmm_           pztrmm
#define pytrmm_           pytrmm
#define pztrsm_           pztrsm
#define pytrsm_           pytrsm
#define pztranu_          pztranu
#define pytranu_          pytranu
#define pztranc_          pztranc
#define pytranc_          pytranc

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )
/*
* These defines set up the naming scheme required to have a FORTRAN
* routine call a C routine (which is what the PBLAS are written in)
* for systems where the fortran "compiler" is actually f2c (a Fortran
* to C conversion utility).
*/

/*
* Initialization routines
*/
#define blacs_pinfo_    blacs_pinfo__
#define blacs_setup_    blacs_setup__
#define blacs_set_      blacs_set__
#define blacs_get_      blacs_get__
#define blacs_gridinit_ blacs_gridinit__
#define blacs_gridmap_  blacs_gridmap__
/*
* Destruction routines
*/
#define blacs_freebuff_ blacs_freebuff__
#define blacs_gridexit_ blacs_gridexit__
#define blacs_abort_    blacs_abort__
#define blacs_exit_     blacs_exit__
/*
* Informational & misc.
*/
#define blacs_gridinfo_ blacs_gridinfo__
#define blacs_pnum_     blacs_pnum__
#define blacs_pcoord_   blacs_pcoord__
#define blacs_barrier_  blacs_barrier__
#endif
