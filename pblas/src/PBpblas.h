/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  This file includes PBLAS definitions. All PBLAS routines include this
*  file.
*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine. No redefinition is necessary  to  have  the
*  following FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALL PDGEMM(...)               void pdgemm_(...)
*           CALL PEGEMM(...)               void pegemm_(...)
*
*  This is the PBLAS default.
*/
#define    PB_freebuf_         PB_freebuf_
#define    PB_topget_          pb_topget_
#define    PB_topset_          pb_topset_

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine. No redefinition is necessary  to  have  the
*  following FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALL PDGEMM(...)               void PDGEMM(...)
*           CALL PEGEMM(...)               void PEGEMM(...)
*/
#define    pilaenv_            PILAENV
#define    PB_freebuf_         PB_FREEBUF
#define    PB_topget_          PB_TOPGET
#define    PB_topset_          PB_TOPSET
                                                     /* Level-1 PBLAS */
#define    picopy_             PICOPY
#define    pscopy_             PSCOPY
#define    pecopy_             PECOPY
#define    pdcopy_             PDCOPY
#define    pecopy_             PECOPY
#define    pccopy_             PCCOPY
#define    pycopy_             PYCOPY
#define    pzcopy_             PZCOPY
#define    pycopy_             PYCOPY

#define    psswap_             PSSWAP
#define    peswap_             PESWAP
#define    pdswap_             PDSWAP
#define    peswap_             PESWAP
#define    pcswap_             PCSWAP
#define    pyswap_             PYSWAP
#define    pzswap_             PZSWAP
#define    pyswap_             PYSWAP

#define    psaxpy_             PSAXPY
#define    peaxpy_             PEAXPY
#define    pdaxpy_             PDAXPY
#define    peaxpy_             PEAXPY
#define    pcaxpy_             PCAXPY
#define    pyaxpy_             PYAXPY
#define    pzaxpy_             PZAXPY
#define    pyaxpy_             PYAXPY

#define    psscal_             PSSCAL
#define    pescal_             PESCAL
#define    pdscal_             PDSCAL
#define    pescal_             PESCAL
#define    pcscal_             PCSCAL
#define    pyscal_             PYSCAL
#define    pzscal_             PZSCAL
#define    pyscal_             PYSCAL
#define    pcsscal_            PCSSCAL
#define    pyescal_            PYESCAL
#define    pzdscal_            PZDSCAL
#define    pyescal_            PYESCAL

#define    psasum_             PSASUM
#define    peasum_             PEASUM
#define    pdasum_             PDASUM
#define    peasum_             PEASUM
#define    pscasum_            PSCASUM
#define    peyasum_            PEYASUM
#define    pdzasum_            PDZASUM
#define    peyasum_            PEYASUM

#define    psnrm2_             PSNRM2
#define    penrm2_             PENRM2
#define    pdnrm2_             PDNRM2
#define    penrm2_             PENRM2
#define    pscnrm2_            PSCNRM2
#define    peynrm2_            PEYNRM2
#define    pdznrm2_            PDZNRM2
#define    peynrm2_            PEYNRM2

#define    psdot_              PSDOT
#define    pedot_              PEDOT
#define    pddot_              PDDOT
#define    pedot_              PEDOT
#define    pcdotu_             PCDOTU
#define    pydotu_             PYDOTU
#define    pzdotu_             PZDOTU
#define    pydotu_             PYDOTU
#define    pcdotc_             PCDOTC
#define    pydotc_             PYDOTC
#define    pzdotc_             PZDOTC
#define    pydotc_             PYDOTC

#define    psamax_             PSAMAX
#define    peamax_             PEAMAX
#define    pdamax_             PDAMAX
#define    peamax_             PEAMAX
#define    pcamax_             PCAMAX
#define    pyamax_             PYAMAX
#define    pzamax_             PZAMAX
#define    pyamax_             PYAMAX

#define    psgemv_             PSGEMV
#define    pegemv_             PEGEMV
#define    pdgemv_             PDGEMV
#define    pegemv_             PEGEMV
#define    pcgemv_             PCGEMV
#define    pygemv_             PYGEMV
#define    pzgemv_             PZGEMV
#define    pygemv_             PYGEMV

#define    psagemv_            PSAGEMV
#define    peagemv_            PEAGEMV
#define    pdagemv_            PDAGEMV
#define    peagemv_            PEAGEMV
#define    pcagemv_            PCAGEMV
#define    pyagemv_            PYAGEMV
#define    pzagemv_            PZAGEMV
#define    pyagemv_            PYAGEMV

#define    pssymv_             PSSYMV
#define    pesymv_             PESYMV
#define    pdsymv_             PDSYMV
#define    pesymv_             PESYMV
#define    pchemv_             PCHEMV
#define    pyhemv_             PYHEMV
#define    pzhemv_             PZHEMV
#define    pyhemv_             PYHEMV

#define    psasymv_            PSASYMV
#define    peasymv_            PEASYMV
#define    pdasymv_            PDASYMV
#define    peasymv_            PEASYMV
#define    pcahemv_            PCAHEMV
#define    pyahemv_            PYAHEMV
#define    pzahemv_            PZAHEMV
#define    pyahemv_            PYAHEMV

#define    pstrmv_             PSTRMV
#define    petrmv_             PETRMV
#define    pdtrmv_             PDTRMV
#define    petrmv_             PETRMV
#define    pctrmv_             PCTRMV
#define    pytrmv_             PYTRMV
#define    pztrmv_             PZTRMV
#define    pytrmv_             PYTRMV

#define    psatrmv_            PSATRMV
#define    peatrmv_            PEATRMV
#define    pdatrmv_            PDATRMV
#define    peatrmv_            PEATRMV
#define    pcatrmv_            PCATRMV
#define    pyatrmv_            PYATRMV
#define    pzatrmv_            PZATRMV
#define    pyatrmv_            PYATRMV

#define    pstrsv_             PSTRSV
#define    petrsv_             PETRSV
#define    pdtrsv_             PDTRSV
#define    petrsv_             PETRSV
#define    pctrsv_             PCTRSV
#define    pytrsv_             PYTRSV
#define    pztrsv_             PZTRSV
#define    pytrsv_             PYTRSV

#define    psger_              PSGER
#define    peger_              PEGER
#define    pdger_              PDGER
#define    peger_              PEGER
#define    pcgeru_             PCGERU
#define    pygeru_             PYGERU
#define    pzgeru_             PZGERU
#define    pygeru_             PYGERU
#define    pcgerc_             PCGERC
#define    pygerc_             PYGERC
#define    pzgerc_             PZGERC
#define    pygerc_             PYGERC

#define    pssyr_              PSSYR
#define    pesyr_              PESYR
#define    pdsyr_              PDSYR
#define    pesyr_              PESYR
#define    pcher_              PCHER
#define    pyher_              PYHER
#define    pzher_              PZHER
#define    pyher_              PYHER

#define    pssyr2_             PSSYR2
#define    pesyr2_             PESYR2
#define    pdsyr2_             PDSYR2
#define    pesyr2_             PESYR2
#define    pcher2_             PCHER2
#define    pyher2_             PYHER2
#define    pzher2_             PZHER2
#define    pyher2_             PYHER2

#define    psgemm_             PSGEMM
#define    pegemm_             PEGEMM
#define    pdgemm_             PDGEMM
#define    pegemm_             PEGEMM
#define    pcgemm_             PCGEMM
#define    pygemm_             PYGEMM
#define    pzgemm_             PZGEMM
#define    pygemm_             PYGEMM

#define    psgeadd_            PSGEADD
#define    pegeadd_            PEGEADD
#define    pdgeadd_            PDGEADD
#define    pegeadd_            PEGEADD
#define    pcgeadd_            PCGEADD
#define    pygeadd_            PYGEADD
#define    pzgeadd_            PZGEADD
#define    pygeadd_            PYGEADD

#define    pssymm_             PSSYMM
#define    pesymm_             PESYMM
#define    pdsymm_             PDSYMM
#define    pesymm_             PESYMM
#define    pcsymm_             PCSYMM
#define    pysymm_             PYSYMM
#define    pchemm_             PCHEMM
#define    pyhemm_             PYHEMM
#define    pzsymm_             PZSYMM
#define    pysymm_             PYSYMM
#define    pzhemm_             PZHEMM
#define    pyhemm_             PYHEMM

#define    pstrmm_             PSTRMM
#define    petrmm_             PETRMM
#define    pdtrmm_             PDTRMM
#define    petrmm_             PETRMM
#define    pctrmm_             PCTRMM
#define    pytrmm_             PYTRMM
#define    pztrmm_             PZTRMM
#define    pytrmm_             PYTRMM

#define    pstrsm_             PSTRSM
#define    petrsm_             PETRSM
#define    pdtrsm_             PDTRSM
#define    petrsm_             PETRSM
#define    pctrsm_             PCTRSM
#define    pytrsm_             PYTRSM
#define    pztrsm_             PZTRSM
#define    pytrsm_             PYTRSM

#define    pssyrk_             PSSYRK
#define    pesyrk_             PESYRK
#define    pdsyrk_             PDSYRK
#define    pesyrk_             PESYRK
#define    pcsyrk_             PCSYRK
#define    pysyrk_             PYSYRK
#define    pcherk_             PCHERK
#define    pyherk_             PYHERK
#define    pzsyrk_             PZSYRK
#define    pysyrk_             PYSYRK
#define    pzherk_             PZHERK
#define    pyherk_             PYHERK

#define    pssyr2k_            PSSYR2K
#define    pesyr2k_            PESYR2K
#define    pdsyr2k_            PDSYR2K
#define    pesyr2k_            PESYR2K
#define    pcsyr2k_            PCSYR2K
#define    pysyr2k_            PYSYR2K
#define    pcher2k_            PCHER2K
#define    pyher2k_            PYHER2K
#define    pzsyr2k_            PZSYR2K
#define    pysyr2k_            PYSYR2K
#define    pzher2k_            PZHER2K
#define    pyher2k_            PYHER2K

#define    pstradd_            PSTRADD
#define    petradd_            PETRADD
#define    pdtradd_            PDTRADD
#define    petradd_            PETRADD
#define    pctradd_            PCTRADD
#define    pytradd_            PYTRADD
#define    pztradd_            PZTRADD
#define    pytradd_            PYTRADD

#define    pstran_             PSTRAN
#define    petran_             PETRAN
#define    pdtran_             PDTRAN
#define    petran_             PETRAN
#define    pctranu_            PCTRANU
#define    pytranu_            PYTRANU
#define    pztranu_            PZTRANU
#define    pytranu_            PYTRANU
#define    pctranc_            PCTRANC
#define    pytranc_            PYTRANC
#define    pztranc_            PZTRANC
#define    pytranc_            PYTRANC

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine call a C routine with the following  FORTRAN to C interface:
*
*           FORTRAN CALL                   C DECLARATION
*           CALLL PDGEMM(...)              void pdgemm(...)
*           CALLL PEGEMM(...)              void pegemm(...)
*/
#define    pilaenv_            pilaenv
#define    PB_freebuf_         PB_freebuf
#define    PB_topget_          pb_topget
#define    PB_topset_          pb_topset

#define    picopy_             picopy
#define    pscopy_             pscopy
#define    pecopy_             pecopy
#define    pdcopy_             pdcopy
#define    pecopy_             pecopy
#define    pccopy_             pccopy
#define    pycopy_             pycopy
#define    pzcopy_             pzcopy
#define    pycopy_             pycopy

#define    psswap_             psswap
#define    peswap_             peswap
#define    pdswap_             pdswap
#define    peswap_             peswap
#define    pcswap_             pcswap
#define    pyswap_             pyswap
#define    pzswap_             pzswap
#define    pyswap_             pyswap

#define    psaxpy_             psaxpy
#define    peaxpy_             peaxpy
#define    pdaxpy_             pdaxpy
#define    peaxpy_             peaxpy
#define    pcaxpy_             pcaxpy
#define    pyaxpy_             pyaxpy
#define    pzaxpy_             pzaxpy
#define    pyaxpy_             pyaxpy

#define    psscal_             psscal
#define    pescal_             pescal
#define    pdscal_             pdscal
#define    pescal_             pescal
#define    pcscal_             pcscal
#define    pyscal_             pyscal
#define    pzscal_             pzscal
#define    pyscal_             pyscal
#define    pcsscal_            pcsscal
#define    pyescal_            pyescal
#define    pzdscal_            pzdscal
#define    pyescal_            pyescal

#define    psasum_             psasum
#define    peasum_             peasum
#define    pdasum_             pdasum
#define    peasum_             peasum
#define    pscasum_            pscasum
#define    peyasum_            peyasum
#define    pdzasum_            pdzasum
#define    peyasum_            peyasum

#define    psnrm2_             psnrm2
#define    penrm2_             penrm2
#define    pdnrm2_             pdnrm2
#define    penrm2_             penrm2
#define    pscnrm2_            pscnrm2
#define    peynrm2_            peynrm2
#define    pdznrm2_            pdznrm2
#define    peynrm2_            peynrm2

#define    psdot_              psdot
#define    pedot_              pedot
#define    pddot_              pddot
#define    pedot_              pedot
#define    pcdotu_             pcdotu
#define    pydotu_             pydotu
#define    pzdotu_             pzdotu
#define    pydotu_             pydotu
#define    pcdotc_             pcdotc
#define    pydotc_             pydotc
#define    pzdotc_             pzdotc
#define    pydotc_             pydotc

#define    psamax_             psamax
#define    peamax_             peamax
#define    pdamax_             pdamax
#define    peamax_             peamax
#define    pcamax_             pcamax
#define    pyamax_             pyamax
#define    pzamax_             pzamax
#define    pyamax_             pyamax

#define    psgemv_             psgemv
#define    pegemv_             pegemv
#define    pdgemv_             pdgemv
#define    pegemv_             pegemv
#define    pcgemv_             pcgemv
#define    pygemv_             pygemv
#define    pzgemv_             pzgemv
#define    pygemv_             pygemv

#define    psagemv_            psagemv
#define    peagemv_            peagemv
#define    pdagemv_            pdagemv
#define    peagemv_            peagemv
#define    pcagemv_            pcagemv
#define    pyagemv_            pyagemv
#define    pzagemv_            pzagemv
#define    pyagemv_            pyagemv

#define    pssymv_             pssymv
#define    pesymv_             pesymv
#define    pdsymv_             pdsymv
#define    pesymv_             pesymv
#define    pchemv_             pchemv
#define    pyhemv_             pyhemv
#define    pzhemv_             pzhemv
#define    pyhemv_             pyhemv

#define    psasymv_            psasymv
#define    peasymv_            peasymv
#define    pdasymv_            pdasymv
#define    peasymv_            peasymv
#define    pcahemv_            pcahemv
#define    pyahemv_            pyahemv
#define    pzahemv_            pzahemv
#define    pyahemv_            pyahemv

#define    pstrmv_             pstrmv
#define    petrmv_             petrmv
#define    pdtrmv_             pdtrmv
#define    petrmv_             petrmv
#define    pctrmv_             pctrmv
#define    pytrmv_             pytrmv
#define    pztrmv_             pztrmv
#define    pytrmv_             pytrmv

#define    psatrmv_            psatrmv
#define    peatrmv_            peatrmv
#define    pdatrmv_            pdatrmv
#define    peatrmv_            peatrmv
#define    pcatrmv_            pcatrmv
#define    pyatrmv_            pyatrmv
#define    pzatrmv_            pzatrmv
#define    pyatrmv_            pyatrmv

#define    pstrsv_             pstrsv
#define    petrsv_             petrsv
#define    pdtrsv_             pdtrsv
#define    petrsv_             petrsv
#define    pctrsv_             pctrsv
#define    pytrsv_             pytrsv
#define    pztrsv_             pztrsv
#define    pytrsv_             pytrsv

#define    psger_              psger
#define    peger_              peger
#define    pdger_              pdger
#define    peger_              peger
#define    pcgeru_             pcgeru
#define    pygeru_             pygeru
#define    pzgeru_             pzgeru
#define    pygeru_             pygeru
#define    pcgerc_             pcgerc
#define    pygerc_             pygerc
#define    pzgerc_             pzgerc
#define    pygerc_             pygerc

#define    pssyr_              pssyr
#define    pesyr_              pesyr
#define    pdsyr_              pdsyr
#define    pesyr_              pesyr
#define    pcher_              pcher
#define    pyher_              pyher
#define    pzher_              pzher
#define    pyher_              pyher

#define    pssyr2_             pssyr2
#define    pesyr2_             pesyr2
#define    pdsyr2_             pdsyr2
#define    pesyr2_             pesyr2
#define    pcher2_             pcher2
#define    pyher2_             pyher2
#define    pzher2_             pzher2
#define    pyher2_             pyher2

#define    psgeadd_            psgeadd
#define    pegeadd_            pegeadd
#define    pdgeadd_            pdgeadd
#define    pegeadd_            pegeadd
#define    pcgeadd_            pcgeadd
#define    pygeadd_            pygeadd
#define    pzgeadd_            pzgeadd
#define    pygeadd_            pygeadd

#define    psgemm_             psgemm
#define    pegemm_             pegemm
#define    pdgemm_             pdgemm
#define    pegemm_             pegemm
#define    pcgemm_             pcgemm
#define    pygemm_             pygemm
#define    pzgemm_             pzgemm
#define    pygemm_             pygemm

#define    pssymm_             pssymm
#define    pesymm_             pesymm
#define    pdsymm_             pdsymm
#define    pesymm_             pesymm
#define    pcsymm_             pcsymm
#define    pysymm_             pysymm
#define    pchemm_             pchemm
#define    pyhemm_             pyhemm
#define    pzsymm_             pzsymm
#define    pysymm_             pysymm
#define    pzhemm_             pzhemm
#define    pyhemm_             pyhemm

#define    pstrmm_             pstrmm
#define    petrmm_             petrmm
#define    pdtrmm_             pdtrmm
#define    petrmm_             petrmm
#define    pctrmm_             pctrmm
#define    pytrmm_             pytrmm
#define    pztrmm_             pztrmm
#define    pytrmm_             pytrmm

#define    pstrsm_             pstrsm
#define    petrsm_             petrsm
#define    pdtrsm_             pdtrsm
#define    petrsm_             petrsm
#define    pctrsm_             pctrsm
#define    pytrsm_             pytrsm
#define    pztrsm_             pztrsm
#define    pytrsm_             pytrsm

#define    pssyrk_             pssyrk
#define    pesyrk_             pesyrk
#define    pdsyrk_             pdsyrk
#define    pesyrk_             pesyrk
#define    pcsyrk_             pcsyrk
#define    pysyrk_             pysyrk
#define    pcherk_             pcherk
#define    pyherk_             pyherk
#define    pzsyrk_             pzsyrk
#define    pysyrk_             pysyrk
#define    pzherk_             pzherk
#define    pyherk_             pyherk

#define    pssyr2k_            pssyr2k
#define    pesyr2k_            pesyr2k
#define    pdsyr2k_            pdsyr2k
#define    pesyr2k_            pesyr2k
#define    pcsyr2k_            pcsyr2k
#define    pysyr2k_            pysyr2k
#define    pcher2k_            pcher2k
#define    pyher2k_            pyher2k
#define    pzsyr2k_            pzsyr2k
#define    pysyr2k_            pysyr2k
#define    pzher2k_            pzher2k
#define    pyher2k_            pyher2k

#define    pstradd_            pstradd
#define    petradd_            petradd
#define    pdtradd_            pdtradd
#define    petradd_            petradd
#define    pctradd_            pctradd
#define    pytradd_            pytradd
#define    pztradd_            pztradd
#define    pytradd_            pytradd

#define    pstran_             pstran
#define    petran_             petran
#define    pdtran_             pdtran
#define    petran_             petran
#define    pctranu_            pctranu
#define    pytranu_            pytranu
#define    pztranu_            pztranu
#define    pytranu_            pytranu
#define    pctranc_            pctranc
#define    pytranc_            pytranc
#define    pztranc_            pztranc
#define    pytranc_            pytranc

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )

#define    PB_freebuf_         PB_freebuf__
#define    PB_topget_          pb_topget__
#define    PB_topset_          pb_topset__

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

void           PB_freebuf_     ( void );

void           PB_topget_      ( Int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           PB_topset_      ( Int *,     F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T );

void           picopy_         ( Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int * );
void           pscopy_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pecopy_         ( Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdcopy_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pecopy_         ( Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );
void           pccopy_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pycopy_         ( Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzcopy_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pycopy_         ( Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );

void           psswap_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           peswap_         ( Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdswap_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           peswap_         ( Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcswap_         ( Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pyswap_         ( Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzswap_         ( Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pyswap_         ( Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );

void           psaxpy_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           peaxpy_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdaxpy_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           peaxpy_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcaxpy_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pyaxpy_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzaxpy_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pyaxpy_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int * );

void           psscal_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pescal_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdscal_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pescal_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcscal_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyscal_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcsscal_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyescal_        ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzscal_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyscal_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzdscal_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyescal_        ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psasum_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peasum_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdasum_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peasum_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pscasum_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peyasum_        ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdzasum_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peyasum_        ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psnrm2_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           penrm2_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdnrm2_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           penrm2_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pscnrm2_        ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peynrm2_        ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdznrm2_        ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peynrm2_        ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psdot_          ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pedot_          ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pddot_          ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pedot_          ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcdotc_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pydotc_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pcdotu_         ( Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pydotu_         ( Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzdotc_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pydotc_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int * );
void           pzdotu_         ( Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pydotu_         ( Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int * );

void           psamax_         ( Int *,     float *,   Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           peamax_         ( Int *,     EREAL *,   Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdamax_         ( Int *,     double *,  Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           peamax_         ( Int *,     EREAL *,  Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcamax_         ( Int *,     float *,   Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pyamax_         ( Int *,     EREAL *,   Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzamax_         ( Int *,     double *,  Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pyamax_         ( Int *,     EREAL *,  Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );

void           psgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pegemv_         ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pegemv_         ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pygemv_         ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzgemv_         ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pygemv_         ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           peagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pcagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pzagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pyagemv_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psger_          ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           peger_          ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int * );
void           pdger_          ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           peger_          ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int * );
void           pcgerc_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pygerc_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int * );
void           pcgeru_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pygeru_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int * );
void           pzgerc_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pygerc_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int * );
void           pzgeru_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pygeru_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int * );

void           pssymv_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pesymv_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdsymv_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pesymv_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     Int * );
void           pchemv_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pyhemv_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzhemv_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pyhemv_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     Int * );

void           psasymv_        ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           peasymv_        ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pdasymv_        ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           peasymv_        ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     Int * );
void           pcahemv_        ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     Int * );
void           pyahemv_        ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     Int * );
void           pzahemv_        ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     Int * );
void           pyahemv_        ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     Int * );

void           pssyr_          ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pesyr_          ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pdsyr_          ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pesyr_          ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pcher_          ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pyher_          ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pzher_          ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pyher_          ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );

void           pssyr2_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pesyr2_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int * );
void           pdsyr2_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pesyr2_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int * );
void           pcher2_         ( F_CHAR_T,  Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int * );
void           pyher2_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int * );
void           pzher2_         ( F_CHAR_T,  Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int * );
void           pyher2_         ( F_CHAR_T,  Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int * );

void           pstrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           petrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdtrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           petrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pctrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pytrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pztrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pytrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           peatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pdatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           peatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );
void           pcatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     Int * );
void           pyatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     Int * );
void           pzatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     Int * );
void           pyatrmv_        ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     Int * );

void           pstrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           petrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pdtrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           petrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pctrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pytrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pztrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           pytrsv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int *,
                                 Int * );

void           psgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pegeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int * );
void           pdgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pegeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int * );
void           pcgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pygeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int * );
void           pzgeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pygeadd_        ( F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int * );

void           psgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pegemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int * );
void           pdgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pegemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int * );
void           pcgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   Int *,
                                 Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int * );
void           pygemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int * );
void           pzgemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  Int *,
                                 Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int * );
void           pygemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int * );

void           pssymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pesymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pdsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pesymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pcsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pysymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pzsymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pysymm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pchemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pyhemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pzhemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pyhemm_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );

void           pssyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pesyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pdsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pesyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pcsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pysyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pzsyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pysyr2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pcher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pyher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pzher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pyher2k_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );

void           pssyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pesyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int * );
void           pdsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pesyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int * );
void           pcsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pysyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int * );
void           pzsyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pysyrk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int * );
void           pcherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pyherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int * );
void           pzherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pyherk_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int * );

void           pstradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           petradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int * );
void           pdtradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           petradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int * );
void           pctradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int * );
void           pytradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int * );
void           pztradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int * );
void           pytradd_        ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int * );

void           pstran_         ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           petran_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pdtran_         ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           petran_         ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pctranc_        ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pytranc_        ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pztranc_        ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pytranc_        ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );
void           pctranu_        ( Int *,     Int *,     float *,
                                 float *,   Int *,     Int *,
                                 Int *,     float *,   float *,
                                 Int *,     Int *,     Int * );
void           pytranu_        ( Int *,     Int *,     EREAL *,
                                 EREAL *,   Int *,     Int *,
                                 Int *,     EREAL *,   EREAL *,
                                 Int *,     Int *,     Int * );
void           pztranu_        ( Int *,     Int *,     double *,
                                 double *,  Int *,     Int *,
                                 Int *,     double *,  double *,
                                 Int *,     Int *,     Int * );
void           pytranu_        ( Int *,     Int *,     EREAL *,
                                 EREAL *,  Int *,     Int *,
                                 Int *,     EREAL *,  EREAL *,
                                 Int *,     Int *,     Int * );

void           pstrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           petrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pdtrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           petrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pctrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pytrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pztrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pytrmm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );

void           pstrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           petrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pdtrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           petrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pctrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 float *,   float *,   Int *,
                                 Int *,     Int *,     float *,
                                 Int *,     Int *,     Int * );
void           pytrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,   EREAL *,   Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
void           pztrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 double *,  double *,  Int *,
                                 Int *,     Int *,     double *,
                                 Int *,     Int *,     Int * );
void           pytrsm_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 EREAL *,  EREAL *,  Int *,
                                 Int *,     Int *,     EREAL *,
                                 Int *,     Int *,     Int * );
#else

void           PB_freebuf_     ();
void           PB_topget_      ();
void           PB_topset_      ();

void           picopy_         ();
void           pscopy_         ();
void           pecopy_         ();
void           pdcopy_         ();
void           pecopy_         ();
void           pccopy_         ();
void           pycopy_         ();
void           pzcopy_         ();
void           pycopy_         ();

void           psswap_         ();
void           peswap_         ();
void           pdswap_         ();
void           peswap_         ();
void           pcswap_         ();
void           pyswap_         ();
void           pzswap_         ();
void           pyswap_         ();

void           psaxpy_         ();
void           peaxpy_         ();
void           pdaxpy_         ();
void           peaxpy_         ();
void           pcaxpy_         ();
void           pyaxpy_         ();
void           pzaxpy_         ();
void           pyaxpy_         ();

void           psscal_         ();
void           pescal_         ();
void           pdscal_         ();
void           pescal_         ();
void           pcscal_         ();
void           pyscal_         ();
void           pcsscal_        ();
void           pyescal_        ();
void           pzscal_         ();
void           pyscal_         ();
void           pzdscal_        ();
void           pyescal_        ();

void           psasum_         ();
void           peasum_         ();
void           pdasum_         ();
void           peasum_         ();
void           pscasum_        ();
void           peyasum_        ();
void           pdzasum_        ();
void           peyasum_        ();

void           psnrm2_         ();
void           penrm2_         ();
void           pdnrm2_         ();
void           penrm2_         ();
void           pscnrm2_        ();
void           peynrm2_        ();
void           pdznrm2_        ();
void           peynrm2_        ();

void           psdot_          ();
void           pedot_          ();
void           pddot_          ();
void           pedot_          ();
void           pcdotc_         ();
void           pydotc_         ();
void           pcdotu_         ();
void           pydotu_         ();
void           pzdotc_         ();
void           pydotc_         ();
void           pzdotu_         ();
void           pydotu_         ();

void           psamax_         ();
void           peamax_         ();
void           pdamax_         ();
void           peamax_         ();
void           pcamax_         ();
void           pyamax_         ();
void           pzamax_         ();
void           pyamax_         ();

void           psgemv_         ();
void           pegemv_         ();
void           pdgemv_         ();
void           pegemv_         ();
void           pcgemv_         ();
void           pygemv_         ();
void           pzgemv_         ();
void           pygemv_         ();

void           psagemv_        ();
void           peagemv_        ();
void           pdagemv_        ();
void           peagemv_        ();
void           pcagemv_        ();
void           pyagemv_        ();
void           pzagemv_        ();
void           pyagemv_        ();

void           psger_          ();
void           peger_          ();
void           pdger_          ();
void           peger_          ();
void           pcgerc_         ();
void           pygerc_         ();
void           pcgeru_         ();
void           pygeru_         ();
void           pzgerc_         ();
void           pygerc_         ();
void           pzgeru_         ();
void           pygeru_         ();

void           pssymv_         ();
void           pesymv_         ();
void           pdsymv_         ();
void           pesymv_         ();
void           pchemv_         ();
void           pyhemv_         ();
void           pzhemv_         ();
void           pyhemv_         ();

void           psasymv_        ();
void           peasymv_        ();
void           pdasymv_        ();
void           peasymv_        ();
void           pcahemv_        ();
void           pyahemv_        ();
void           pzahemv_        ();
void           pyahemv_        ();

void           pssyr_          ();
void           pesyr_          ();
void           pdsyr_          ();
void           pesyr_          ();
void           pcher_          ();
void           pyher_          ();
void           pzher_          ();
void           pyher_          ();

void           pssyr2_         ();
void           pesyr2_         ();
void           pdsyr2_         ();
void           pesyr2_         ();
void           pcher2_         ();
void           pyher2_         ();
void           pzher2_         ();
void           pyher2_         ();

void           pstrmv_         ();
void           petrmv_         ();
void           pdtrmv_         ();
void           petrmv_         ();
void           pctrmv_         ();
void           pytrmv_         ();
void           pztrmv_         ();
void           pytrmv_         ();

void           psatrmv_        ();
void           peatrmv_        ();
void           pdatrmv_        ();
void           peatrmv_        ();
void           pcatrmv_        ();
void           pyatrmv_        ();
void           pzatrmv_        ();
void           pyatrmv_        ();

void           pstrsv_         ();
void           petrsv_         ();
void           pdtrsv_         ();
void           petrsv_         ();
void           pctrsv_         ();
void           pytrsv_         ();
void           pztrsv_         ();
void           pytrsv_         ();

void           psgeadd_        ();
void           pegeadd_        ();
void           pdgeadd_        ();
void           pegeadd_        ();
void           pcgeadd_        ();
void           pygeadd_        ();
void           pzgeadd_        ();
void           pygeadd_        ();

void           psgemm_         ();
void           pegemm_         ();
void           pdgemm_         ();
void           pegemm_         ();
void           pcgemm_         ();
void           pygemm_         ();
void           pzgemm_         ();
void           pygemm_         ();

void           pssymm_         ();
void           pesymm_         ();
void           pdsymm_         ();
void           pesymm_         ();
void           pcsymm_         ();
void           pysymm_         ();
void           pchemm_         ();
void           pyhemm_         ();
void           pzsymm_         ();
void           pysymm_         ();
void           pzhemm_         ();
void           pyhemm_         ();

void           pssyr2k_        ();
void           pesyr2k_        ();
void           pdsyr2k_        ();
void           pesyr2k_        ();
void           pcsyr2k_        ();
void           pysyr2k_        ();
void           pcher2k_        ();
void           pyher2k_        ();
void           pzsyr2k_        ();
void           pysyr2k_        ();
void           pzher2k_        ();
void           pyher2k_        ();

void           pssyrk_         ();
void           pesyrk_         ();
void           pdsyrk_         ();
void           pesyrk_         ();
void           pcsyrk_         ();
void           pysyrk_         ();
void           pcherk_         ();
void           pyherk_         ();
void           pzsyrk_         ();
void           pysyrk_         ();
void           pzherk_         ();
void           pyherk_         ();

void           pstradd_        ();
void           petradd_        ();
void           pdtradd_        ();
void           petradd_        ();
void           pctradd_        ();
void           pytradd_        ();
void           pztradd_        ();
void           pytradd_        ();

void           pstran_         ();
void           petran_         ();
void           pdtran_         ();
void           petran_         ();
void           pctranc_        ();
void           pytranc_        ();
void           pctranu_        ();
void           pytranu_        ();
void           pztranc_        ();
void           pytranc_        ();
void           pztranu_        ();
void           pytranu_        ();

void           pstrmm_         ();
void           petrmm_         ();
void           pdtrmm_         ();
void           petrmm_         ();
void           pctrmm_         ();
void           pytrmm_         ();
void           pztrmm_         ();
void           pytrmm_         ();

void           pstrsm_         ();
void           petrsm_         ();
void           pdtrsm_         ();
void           petrsm_         ();
void           pctrsm_         ();
void           pytrsm_         ();
void           pztrsm_         ();
void           pytrsm_         ();

#endif
