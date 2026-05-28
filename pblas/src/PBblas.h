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
*  This file includes F77 BLAS definitions. All  PBLAS  routines include
*  this file.
*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
#define    CNOTRAN             'N'
#define    CNOCONJG            'N'
#define    CTRAN               'T'
#define    CCONJG              'Z'
#define    CCOTRAN             'C'

#define    CALL                'A'
#define    CLOWER              'L'
#define    CUPPER              'U'
#define    CDIAGONAL           'D'

#define    CLEFT               'L'
#define    CRIGHT              'R'

#define    CUNIT               'U'
#define    CNOUNIT             'N'

#define    CINIT               'I'
#define    CNOINIT             'N'

#define    CFORWARD            'F'
#define    CBACKWARD           'B'

#define    CREUSE              'R'
#define    CALLOCATE           'A'

#define    NOTRAN              "N"
#define    NOCONJG             "N"
#define    TRAN                "T"
#define    CONJG               "Z"
#define    COTRAN              "C"

#define    ALL                 "A"
#define    LOWER               "L"
#define    UPPER               "U"
#define    DIAGONAL            "D"

#define    LEFT                "L"
#define    RIGHT               "R"

#define    UNIT                "U"
#define    NOUNIT              "N"

#define    INIT                "I"
#define    NOINIT              "N"

#define    FORWARD             "F"
#define    BACKWARD            "B"

#define    REUSE               "R"
#define    ALLOCATE            "A"

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine. No redefinition is necessary  to  have
*  the following FORTRAN to C interface:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          dgemm_(...)
*           SUBROUTINE EGEMM(...)          egemm_(...)
*
*  This is the PBLAS default.
*/
#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          DGEMM(...)
*           SUBROUTINE EGEMM(...)          EGEMM(...)
*/
#define    srot_               SROT
#define    erot_               EROT
#define    drot_               DROT
#define    erot_               EROT

#define    sswap_              SSWAP
#define    eswap_              ESWAP
#define    dswap_              DSWAP
#define    eswap_              ESWAP
#define    cswap_              CSWAP
#define    yswap_              YSWAP
#define    zswap_              ZSWAP
#define    yswap_              YSWAP

#define    scopy_              SCOPY
#define    ecopy_              ECOPY
#define    dcopy_              DCOPY
#define    ecopy_              ECOPY
#define    ccopy_              CCOPY
#define    ycopy_              YCOPY
#define    zcopy_              ZCOPY
#define    ycopy_              YCOPY

#define    saxpy_              SAXPY
#define    eaxpy_              EAXPY
#define    daxpy_              DAXPY
#define    eaxpy_              EAXPY
#define    caxpy_              CAXPY
#define    yaxpy_              YAXPY
#define    zaxpy_              ZAXPY
#define    yaxpy_              YAXPY

#define    sscal_              SSCAL
#define    escal_              ESCAL
#define    dscal_              DSCAL
#define    escal_              ESCAL
#define    cscal_              CSCAL
#define    yscal_              YSCAL
#define    zscal_              ZSCAL
#define    yscal_              YSCAL
#define    csscal_             CSSCAL
#define    yescal_             YESCAL
#define    zdscal_             ZDSCAL
#define    yescal_             YESCAL

#define    sasum_              SASUM
#define    easum_              EASUM
#define    dasum_              DASUM
#define    easum_              EASUM
#define    scasum_             SCASUM
#define    eyasum_             EYASUM
#define    dzasum_             DZASUM
#define    eyasum_             EYASUM

#define    snrm2_              SNRM2
#define    enrm2_              ENRM2
#define    dnrm2_              DNRM2
#define    enrm2_              ENRM2
#define    scnrm2_             SCNRM2
#define    eynrm2_             EYNRM2
#define    dznrm2_             DZNRM2
#define    eynrm2_             EYNRM2

#define    sdot_               SDOT
#define    edot_               EDOT
#define    ddot_               DDOT
#define    edot_               EDOT
#define    cdotu_              CDOTU
#define    ydotu_              YDOTU
#define    zdotu_              ZDOTU
#define    ydotu_              YDOTU
#define    cdotc_              CDOTC
#define    ydotc_              YDOTC
#define    zdotc_              ZDOTC
#define    ydotc_              YDOTC

#define    isamax_             ISAMAX
#define    ieamax_             IEAMAX
#define    idamax_             IDAMAX
#define    ieamax_             IEAMAX
#define    icamax_             ICAMAX
#define    iyamax_             IYAMAX
#define    izamax_             IZAMAX
#define    iyamax_             IYAMAX

#define    sgemv_              SGEMV
#define    egemv_              EGEMV
#define    dgemv_              DGEMV
#define    egemv_              EGEMV
#define    cgemv_              CGEMV
#define    ygemv_              YGEMV
#define    zgemv_              ZGEMV
#define    ygemv_              YGEMV

#define    ssymv_              SSYMV
#define    esymv_              ESYMV
#define    dsymv_              DSYMV
#define    esymv_              ESYMV
#define    chemv_              CHEMV
#define    yhemv_              YHEMV
#define    zhemv_              ZHEMV
#define    yhemv_              YHEMV

#define    strmv_              STRMV
#define    etrmv_              ETRMV
#define    dtrmv_              DTRMV
#define    etrmv_              ETRMV
#define    ctrmv_              CTRMV
#define    ytrmv_              YTRMV
#define    ztrmv_              ZTRMV
#define    ytrmv_              YTRMV

#define    strsv_              STRSV
#define    etrsv_              ETRSV
#define    dtrsv_              DTRSV
#define    etrsv_              ETRSV
#define    ctrsv_              CTRSV
#define    ytrsv_              YTRSV
#define    ztrsv_              ZTRSV
#define    ytrsv_              YTRSV

#define    sger_               SGER
#define    eger_               EGER
#define    dger_               DGER
#define    eger_               EGER
#define    cgeru_              CGERU
#define    ygeru_              YGERU
#define    zgeru_              ZGERU
#define    ygeru_              YGERU
#define    cgerc_              CGERC
#define    ygerc_              YGERC
#define    zgerc_              ZGERC
#define    ygerc_              YGERC

#define    ssyr_               SSYR
#define    esyr_               ESYR
#define    dsyr_               DSYR
#define    esyr_               ESYR
#define    cher_               CHER
#define    yher_               YHER
#define    zher_               ZHER
#define    yher_               YHER

#define    ssyr2_              SSYR2
#define    esyr2_              ESYR2
#define    dsyr2_              DSYR2
#define    esyr2_              ESYR2
#define    cher2_              CHER2
#define    yher2_              YHER2
#define    zher2_              ZHER2
#define    yher2_              YHER2

#define    sgemm_              SGEMM
#define    egemm_              EGEMM
#define    dgemm_              DGEMM
#define    egemm_              EGEMM
#define    cgemm_              CGEMM
#define    ygemm_              YGEMM
#define    zgemm_              ZGEMM
#define    ygemm_              YGEMM

#define    ssymm_              SSYMM
#define    esymm_              ESYMM
#define    dsymm_              DSYMM
#define    esymm_              ESYMM
#define    csymm_              CSYMM
#define    ysymm_              YSYMM
#define    chemm_              CHEMM
#define    yhemm_              YHEMM
#define    zsymm_              ZSYMM
#define    ysymm_              YSYMM
#define    zhemm_              ZHEMM
#define    yhemm_              YHEMM

#define    strmm_              STRMM
#define    etrmm_              ETRMM
#define    dtrmm_              DTRMM
#define    etrmm_              ETRMM
#define    ctrmm_              CTRMM
#define    ytrmm_              YTRMM
#define    ztrmm_              ZTRMM
#define    ytrmm_              YTRMM

#define    strsm_              STRSM
#define    etrsm_              ETRSM
#define    dtrsm_              DTRSM
#define    etrsm_              ETRSM
#define    ctrsm_              CTRSM
#define    ytrsm_              YTRSM
#define    ztrsm_              ZTRSM
#define    ytrsm_              YTRSM

#define    ssyrk_              SSYRK
#define    esyrk_              ESYRK
#define    dsyrk_              DSYRK
#define    esyrk_              ESYRK
#define    csyrk_              CSYRK
#define    ysyrk_              YSYRK
#define    cherk_              CHERK
#define    yherk_              YHERK
#define    zsyrk_              ZSYRK
#define    ysyrk_              YSYRK
#define    zherk_              ZHERK
#define    yherk_              YHERK

#define    ssyr2k_             SSYR2K
#define    esyr2k_             ESYR2K
#define    dsyr2k_             DSYR2K
#define    esyr2k_             ESYR2K
#define    csyr2k_             CSYR2K
#define    ysyr2k_             YSYR2K
#define    cher2k_             CHER2K
#define    yher2k_             YHER2K
#define    zsyr2k_             ZSYR2K
#define    ysyr2k_             YSYR2K
#define    zher2k_             ZHER2K
#define    yher2k_             YHER2K

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE DGEMM(...)          dgemm(...)
*           SUBROUTINE EGEMM(...)          egemm(...)
*/
#define    srot_               srot
#define    erot_               erot
#define    drot_               drot
#define    erot_               erot

#define    sswap_              sswap
#define    eswap_              eswap
#define    dswap_              dswap
#define    eswap_              eswap
#define    cswap_              cswap
#define    yswap_              yswap
#define    zswap_              zswap
#define    yswap_              yswap

#define    scopy_              scopy
#define    ecopy_              ecopy
#define    dcopy_              dcopy
#define    ecopy_              ecopy
#define    ccopy_              ccopy
#define    ycopy_              ycopy
#define    zcopy_              zcopy
#define    ycopy_              ycopy

#define    saxpy_              saxpy
#define    eaxpy_              eaxpy
#define    daxpy_              daxpy
#define    eaxpy_              eaxpy
#define    caxpy_              caxpy
#define    yaxpy_              yaxpy
#define    zaxpy_              zaxpy
#define    yaxpy_              yaxpy

#define    sscal_              sscal
#define    escal_              escal
#define    dscal_              dscal
#define    escal_              escal
#define    cscal_              cscal
#define    yscal_              yscal
#define    zscal_              zscal
#define    yscal_              yscal
#define    csscal_             csscal
#define    yescal_             yescal
#define    zdscal_             zdscal
#define    yescal_             yescal

#define    sasum_              sasum
#define    easum_              easum
#define    dasum_              dasum
#define    easum_              easum
#define    scasum_             scasum
#define    eyasum_             eyasum
#define    dzasum_             dzasum
#define    eyasum_             eyasum

#define    snrm2_              snrm2
#define    enrm2_              enrm2
#define    dnrm2_              dnrm2
#define    enrm2_              enrm2
#define    scnrm2_             scnrm2
#define    eynrm2_             eynrm2
#define    dznrm2_             dznrm2
#define    eynrm2_             eynrm2

#define    sdot_               sdot
#define    edot_               edot
#define    ddot_               ddot
#define    edot_               edot
#define    cdotu_              cdotu
#define    ydotu_              ydotu
#define    zdotu_              zdotu
#define    ydotu_              ydotu
#define    cdotc_              cdotc
#define    ydotc_              ydotc
#define    zdotc_              zdotc
#define    ydotc_              ydotc

#define    isamax_             isamax
#define    ieamax_             ieamax
#define    idamax_             idamax
#define    ieamax_             ieamax
#define    icamax_             icamax
#define    iyamax_             iyamax
#define    izamax_             izamax
#define    iyamax_             iyamax

#define    sgemv_              sgemv
#define    egemv_              egemv
#define    dgemv_              dgemv
#define    egemv_              egemv
#define    cgemv_              cgemv
#define    ygemv_              ygemv
#define    zgemv_              zgemv
#define    ygemv_              ygemv

#define    ssymv_              ssymv
#define    esymv_              esymv
#define    dsymv_              dsymv
#define    esymv_              esymv
#define    chemv_              chemv
#define    yhemv_              yhemv
#define    zhemv_              zhemv
#define    yhemv_              yhemv

#define    strmv_              strmv
#define    etrmv_              etrmv
#define    dtrmv_              dtrmv
#define    etrmv_              etrmv
#define    ctrmv_              ctrmv
#define    ytrmv_              ytrmv
#define    ztrmv_              ztrmv
#define    ytrmv_              ytrmv

#define    strsv_              strsv
#define    etrsv_              etrsv
#define    dtrsv_              dtrsv
#define    etrsv_              etrsv
#define    ctrsv_              ctrsv
#define    ytrsv_              ytrsv
#define    ztrsv_              ztrsv
#define    ytrsv_              ytrsv

#define    sger_               sger
#define    eger_               eger
#define    dger_               dger
#define    eger_               eger
#define    cgeru_              cgeru
#define    ygeru_              ygeru
#define    zgeru_              zgeru
#define    ygeru_              ygeru
#define    cgerc_              cgerc
#define    ygerc_              ygerc
#define    zgerc_              zgerc
#define    ygerc_              ygerc

#define    ssyr_               ssyr
#define    esyr_               esyr
#define    dsyr_               dsyr
#define    esyr_               esyr
#define    cher_               cher
#define    yher_               yher
#define    zher_               zher
#define    yher_               yher

#define    ssyr2_              ssyr2
#define    esyr2_              esyr2
#define    dsyr2_              dsyr2
#define    esyr2_              esyr2
#define    cher2_              cher2
#define    yher2_              yher2
#define    zher2_              zher2
#define    yher2_              yher2

#define    sgemm_              sgemm
#define    egemm_              egemm
#define    dgemm_              dgemm
#define    egemm_              egemm
#define    cgemm_              cgemm
#define    ygemm_              ygemm
#define    zgemm_              zgemm
#define    ygemm_              ygemm

#define    ssymm_              ssymm
#define    esymm_              esymm
#define    dsymm_              dsymm
#define    esymm_              esymm
#define    csymm_              csymm
#define    ysymm_              ysymm
#define    chemm_              chemm
#define    yhemm_              yhemm
#define    zsymm_              zsymm
#define    ysymm_              ysymm
#define    zhemm_              zhemm
#define    yhemm_              yhemm

#define    strmm_              strmm
#define    etrmm_              etrmm
#define    dtrmm_              dtrmm
#define    etrmm_              etrmm
#define    ctrmm_              ctrmm
#define    ytrmm_              ytrmm
#define    ztrmm_              ztrmm
#define    ytrmm_              ytrmm

#define    strsm_              strsm
#define    etrsm_              etrsm
#define    dtrsm_              dtrsm
#define    etrsm_              etrsm
#define    ctrsm_              ctrsm
#define    ytrsm_              ytrsm
#define    ztrsm_              ztrsm
#define    ytrsm_              ytrsm

#define    ssyrk_              ssyrk
#define    esyrk_              esyrk
#define    dsyrk_              dsyrk
#define    esyrk_              esyrk
#define    csyrk_              csyrk
#define    ysyrk_              ysyrk
#define    cherk_              cherk
#define    yherk_              yherk
#define    zsyrk_              zsyrk
#define    ysyrk_              ysyrk
#define    zherk_              zherk
#define    yherk_              yherk

#define    ssyr2k_             ssyr2k
#define    esyr2k_             esyr2k
#define    dsyr2k_             dsyr2k
#define    esyr2k_             esyr2k
#define    csyr2k_             csyr2k
#define    ysyr2k_             ysyr2k
#define    cher2k_             cher2k
#define    yher2k_             yher2k
#define    zsyr2k_             zsyr2k
#define    ysyr2k_             ysyr2k
#define    zher2k_             zher2k
#define    yher2k_             yher2k

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

Int            isamax_         ( Int *,     char *,    Int * );
Int            ieamax_         ( Int *,     char *,    Int * );
Int            idamax_         ( Int *,     char *,    Int * );
Int            ieamax_         ( Int *,     char *,    Int * );
Int            icamax_         ( Int *,     char *,    Int * );
Int            iyamax_         ( Int *,     char *,    Int * );
Int            izamax_         ( Int *,     char *,    Int * );
Int            iyamax_         ( Int *,     char *,    Int * );

F_VOID_FCT     saxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     eaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     daxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     eaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     caxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yaxpy_          ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     scopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ecopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dcopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ecopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ccopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ycopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     zcopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ycopy_          ( Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     escal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     escal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     csscal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yescal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zdscal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yescal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zscal_          ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yscal_          ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     eswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     eswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     cswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     yswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     zswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     yswap_          ( Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     egemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     egemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ygemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zgemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ygemv_          ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     ssymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     esymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dsymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     esymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     chemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yhemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zhemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yhemv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     strmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrmv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     strsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrsv_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     sger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     eger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     eger_           ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cgerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ygerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cgeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ygeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zgerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ygerc_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zgeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ygeru_          ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     ssyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     esyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     dsyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     esyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     cher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     yher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     zher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     yher_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );

F_VOID_FCT     ssyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     esyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dsyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     esyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yher2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     sgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     egemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     egemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ygemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zgemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ygemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ssymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     esymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dsymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     esymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     csymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ysymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zsymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ysymm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     chemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yhemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zhemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yhemm_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     ssyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     esyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dsyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     esyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     csyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ysyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zsyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ysyrk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yherk_          ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ssyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     esyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dsyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     esyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     csyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ysyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zsyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ysyr2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     cher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yher2k_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     strmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrmm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );

F_VOID_FCT     strsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     dtrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     etrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ctrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ztrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );
F_VOID_FCT     ytrsm_          ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int * );

#else

Int            isamax_         ();
Int            ieamax_         ();
Int            idamax_         ();
Int            ieamax_         ();
Int            icamax_         ();
Int            iyamax_         ();
Int            izamax_         ();
Int            iyamax_         ();

F_VOID_FCT     saxpy_          ();
F_VOID_FCT     eaxpy_          ();
F_VOID_FCT     daxpy_          ();
F_VOID_FCT     eaxpy_          ();
F_VOID_FCT     caxpy_          ();
F_VOID_FCT     yaxpy_          ();
F_VOID_FCT     zaxpy_          ();
F_VOID_FCT     yaxpy_          ();

F_VOID_FCT     scopy_          ();
F_VOID_FCT     ecopy_          ();
F_VOID_FCT     dcopy_          ();
F_VOID_FCT     ecopy_          ();
F_VOID_FCT     ccopy_          ();
F_VOID_FCT     ycopy_          ();
F_VOID_FCT     zcopy_          ();
F_VOID_FCT     ycopy_          ();

F_VOID_FCT     sscal_          ();
F_VOID_FCT     escal_          ();
F_VOID_FCT     dscal_          ();
F_VOID_FCT     escal_          ();
F_VOID_FCT     cscal_          ();
F_VOID_FCT     yscal_          ();
F_VOID_FCT     csscal_         ();
F_VOID_FCT     yescal_         ();
F_VOID_FCT     zscal_          ();
F_VOID_FCT     yscal_          ();
F_VOID_FCT     zdscal_         ();
F_VOID_FCT     yescal_         ();

F_VOID_FCT     sswap_          ();
F_VOID_FCT     eswap_          ();
F_VOID_FCT     dswap_          ();
F_VOID_FCT     eswap_          ();
F_VOID_FCT     cswap_          ();
F_VOID_FCT     yswap_          ();
F_VOID_FCT     zswap_          ();
F_VOID_FCT     yswap_          ();

F_VOID_FCT     sgemv_          ();
F_VOID_FCT     egemv_          ();
F_VOID_FCT     dgemv_          ();
F_VOID_FCT     egemv_          ();
F_VOID_FCT     cgemv_          ();
F_VOID_FCT     ygemv_          ();
F_VOID_FCT     zgemv_          ();
F_VOID_FCT     ygemv_          ();

F_VOID_FCT     ssymv_          ();
F_VOID_FCT     esymv_          ();
F_VOID_FCT     dsymv_          ();
F_VOID_FCT     esymv_          ();
F_VOID_FCT     chemv_          ();
F_VOID_FCT     yhemv_          ();
F_VOID_FCT     zhemv_          ();
F_VOID_FCT     yhemv_          ();

F_VOID_FCT     strmv_          ();
F_VOID_FCT     etrmv_          ();
F_VOID_FCT     dtrmv_          ();
F_VOID_FCT     etrmv_          ();
F_VOID_FCT     ctrmv_          ();
F_VOID_FCT     ytrmv_          ();
F_VOID_FCT     ztrmv_          ();
F_VOID_FCT     ytrmv_          ();

F_VOID_FCT     strsv_          ();
F_VOID_FCT     etrsv_          ();
F_VOID_FCT     dtrsv_          ();
F_VOID_FCT     etrsv_          ();
F_VOID_FCT     ctrsv_          ();
F_VOID_FCT     ytrsv_          ();
F_VOID_FCT     ztrsv_          ();
F_VOID_FCT     ytrsv_          ();

F_VOID_FCT     sger_           ();
F_VOID_FCT     eger_           ();
F_VOID_FCT     dger_           ();
F_VOID_FCT     eger_           ();
F_VOID_FCT     cgerc_          ();
F_VOID_FCT     ygerc_          ();
F_VOID_FCT     cgeru_          ();
F_VOID_FCT     ygeru_          ();
F_VOID_FCT     zgerc_          ();
F_VOID_FCT     ygerc_          ();
F_VOID_FCT     zgeru_          ();
F_VOID_FCT     ygeru_          ();

F_VOID_FCT     ssyr_           ();
F_VOID_FCT     esyr_           ();
F_VOID_FCT     dsyr_           ();
F_VOID_FCT     esyr_           ();
F_VOID_FCT     cher_           ();
F_VOID_FCT     yher_           ();
F_VOID_FCT     zher_           ();
F_VOID_FCT     yher_           ();

F_VOID_FCT     ssyr2_          ();
F_VOID_FCT     esyr2_          ();
F_VOID_FCT     dsyr2_          ();
F_VOID_FCT     esyr2_          ();
F_VOID_FCT     cher2_          ();
F_VOID_FCT     yher2_          ();
F_VOID_FCT     zher2_          ();
F_VOID_FCT     yher2_          ();

F_VOID_FCT     sgemm_          ();
F_VOID_FCT     egemm_          ();
F_VOID_FCT     dgemm_          ();
F_VOID_FCT     egemm_          ();
F_VOID_FCT     cgemm_          ();
F_VOID_FCT     ygemm_          ();
F_VOID_FCT     zgemm_          ();
F_VOID_FCT     ygemm_          ();

F_VOID_FCT     ssymm_          ();
F_VOID_FCT     esymm_          ();
F_VOID_FCT     dsymm_          ();
F_VOID_FCT     esymm_          ();
F_VOID_FCT     csymm_          ();
F_VOID_FCT     ysymm_          ();
F_VOID_FCT     zsymm_          ();
F_VOID_FCT     ysymm_          ();
F_VOID_FCT     chemm_          ();
F_VOID_FCT     yhemm_          ();
F_VOID_FCT     zhemm_          ();
F_VOID_FCT     yhemm_          ();

F_VOID_FCT     ssyrk_          ();
F_VOID_FCT     esyrk_          ();
F_VOID_FCT     dsyrk_          ();
F_VOID_FCT     esyrk_          ();
F_VOID_FCT     csyrk_          ();
F_VOID_FCT     ysyrk_          ();
F_VOID_FCT     zsyrk_          ();
F_VOID_FCT     ysyrk_          ();
F_VOID_FCT     cherk_          ();
F_VOID_FCT     yherk_          ();
F_VOID_FCT     zherk_          ();
F_VOID_FCT     yherk_          ();

F_VOID_FCT     ssyr2k_         ();
F_VOID_FCT     esyr2k_         ();
F_VOID_FCT     dsyr2k_         ();
F_VOID_FCT     esyr2k_         ();
F_VOID_FCT     csyr2k_         ();
F_VOID_FCT     ysyr2k_         ();
F_VOID_FCT     zsyr2k_         ();
F_VOID_FCT     ysyr2k_         ();
F_VOID_FCT     cher2k_         ();
F_VOID_FCT     yher2k_         ();
F_VOID_FCT     zher2k_         ();
F_VOID_FCT     yher2k_         ();

F_VOID_FCT     strmm_          ();
F_VOID_FCT     etrmm_          ();
F_VOID_FCT     dtrmm_          ();
F_VOID_FCT     etrmm_          ();
F_VOID_FCT     ctrmm_          ();
F_VOID_FCT     ytrmm_          ();
F_VOID_FCT     ztrmm_          ();
F_VOID_FCT     ytrmm_          ();

F_VOID_FCT     strsm_          ();
F_VOID_FCT     etrsm_          ();
F_VOID_FCT     dtrsm_          ();
F_VOID_FCT     etrsm_          ();
F_VOID_FCT     ctrsm_          ();
F_VOID_FCT     ytrsm_          ();
F_VOID_FCT     ztrsm_          ();
F_VOID_FCT     ytrsm_          ();

#endif
