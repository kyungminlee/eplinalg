/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     March 12, 2002 
*
*  ---------------------------------------------------------------------
*/
/*
*  This file includes PBLAS tools definitions. All PBLAS routines include
*  this file.
*
* ----------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*
*  Descriptor entries for type 1
*/
#define    BLOCK_CYCLIC_2D     1

#define    DTYPE1_             0                   /* Descriptor Type */
#define    CTXT1_              1                     /* BLACS context */
#define    M1_                 2             /* Global Number of Rows */
#define    N1_                 3          /* Global Number of Columns */
#define    MB1_                4                 /* Row Blocking Size */
#define    NB1_                5              /* Column Blocking Size */
#define    RSRC1_              6            /* Starting Processor Row */
#define    CSRC1_              7         /* Starting Processor Column */
#define    LLD1_               8           /* Local Leading Dimension */
#define    DLEN1_              9                 /* Descriptor Length */
/*
*  Descriptor entries for type 2
*/
#define    BLOCK_CYCLIC_2D_INB 2

#define    DTYPE_              0                   /* Descriptor Type */
#define    CTXT_               1                     /* BLACS context */
#define    M_                  2             /* Global Number of Rows */
#define    N_                  3          /* Global Number of Columns */
#define    IMB_                4         /* Initial Row Blocking Size */
#define    INB_                5      /* Initial Column Blocking Size */
#define    MB_                 6                 /* Row Blocking Size */
#define    NB_                 7              /* Column Blocking Size */
#define    RSRC_               8              /* Starting Process Row */
#define    CSRC_               9           /* Starting Process Column */
#define    LLD_                10          /* Local Leading Dimension */
#define    DLEN_               11                /* Descriptor Length */

#define    CPACKING            'P'
#define    CUNPACKING          'U'

#define    PACKING             "P"
#define    UNPACKING           "U"

#define    CGENERAL            'G'
/* #define    CSYMM               'S'  */
/* #define    YSYMM               'S'  */
#define    CHERM               'H'

#define    GENERAL             "G"
#define    SYMM                "S"
#define    HERM                "H"

#define    ONE                 1.0
#define    TWO                 2.0
#define    ZERO                0.0
                            /* Input error checking related constants */
#define    DESCMULT            100
#define    BIGNUM              10000
/*
*  ---------------------------------------------------------------------
*  #define macro functions
*  ---------------------------------------------------------------------
*/
#define    ABS( a_ )           ( ( (a_) <   0  ) ? -(a_) : (a_) )
#define    MIN( a_, b_ )       ( ( (a_) < (b_) ) ?  (a_) : (b_) )
#define    MAX( a_, b_ )       ( ( (a_) > (b_) ) ?  (a_) : (b_) )

#define    FLOOR(a,b) (((a)>0) ? (((a)/(b))) : (-(((-(a))+(b)-1)/(b))))
#define    CEIL(a,b)           ( ( (a)+(b)-1 ) / (b) )
#define    ICEIL(a,b) (((a)>0) ? ((((a)+(b)-1)/(b))) : (-((-(a))/(b))))

#define    Mupcase(C)          (((C)>96 && (C)<123) ? (C) & 0xDF : (C))
#define    Mlowcase(C)         (((C)>64 && (C)< 91) ? (C) | 32   : (C))
/*
*  The following macros perform common modulo operations;  All functions
*  except MPosMod assume arguments are < d (i.e., arguments are themsel-
*  ves within modulo range).
*/
                                                /* increment with mod */
#define    MModInc(I, d)       if(++(I) == (d)) (I) = 0
                                                /* decrement with mod */
#define    MModDec(I, d)       if(--(I) == -1) (I) = (d)-1
                                                   /* positive modulo */
#define    MPosMod(I, d)       ( (I) - ((I)/(d))*(d) )
                                                   /* add two numbers */
#define    MModAdd(I1, I2, d) \
           ( ( (I1) + (I2) < (d) ) ? (I1) + (I2) : (I1) + (I2) - (d) )
                                                        /* add 1 to # */
#define    MModAdd1(I, d) ( ((I) != (d)-1) ? (I) + 1 : 0 )
                                              /* subtract two numbers */
#define    MModSub(I1, I2, d) \
           ( ( (I1) < (I2) ) ? (d) + (I1) - (I2) : (I1) - (I2) )
                                                      /* sub 1 from # */
#define    MModSub1(I, d) ( ((I)!=0) ? (I)-1 : (d)-1 )
/*
*  DNROC computes maximum number of local rows or columns. This macro is
*  only used to compute the time estimates in the Level 3 PBLAS routines.
*/

#define    DNROC( n_, nb_, p_ ) \
           ((double)(((((n_)+(nb_)-1)/(nb_))+(p_)-1)/(p_))*(double)((nb_)))
/*
*  Mptr returns a pointer to a_( i_, j_ ) for readability reasons and
*  also less silly errors ...
*
*  There was some problems with the previous code which read:
*
*      #define    Mptr( a_, i_, j_, lda_, siz_ ) \
*                    ( (a_) + ( ( (i_)+(j_)*(lda_) )*(siz_) ) )
* 
*  since it can overflow the 32-bit integer "easily".
*  The following code should fix the problem.
*  It uses the "off_t" command.
*
*  Change made by Julien Langou on Sat. September 12, 2009. 
*  Fix provided by John Moyard from CNES.
*
*  JL :April 2011: Change off_t by long long
*  off_t is not supported under Windows
*/
#define    Mptr( a_, i_, j_, lda_, siz_ ) \
              ( (a_) + ( (long long) ( (long long)(i_)+ \
              (long long)(j_)*(long long)(lda_))*(long long)(siz_) ) )
/*
*  Mfirstnb and Mlastnb compute the global size of the first and last
*  block corresponding to the interval i_:i_+n_-1 of global indexes.
*/
#define    Mfirstnb( inbt_, n_, i_, inb_, nb_ ) \
              inbt_ = (inb_) - (i_); \
              if( inbt_ <= 0 ) \
                 inbt_ = ( (-inbt_) / (nb_) + 1 ) * (nb_) + inbt_; \
              inbt_ = MIN( inbt_, (n_) )

#define    Mlastnb( inbt_, n_, i_, inb_, nb_ ) \
              inbt_ = (i_) + (n_) - (inb_); \
              if( inbt_ > 0 ) \
              { \
                 inbt_ = -( ( (nb_)+inbt_-1 )/(nb_)-1 )*(nb_) + inbt_; \
                 inbt_ = MIN( inbt_, (n_) ); \
              } \
              else { inbt_ = (n_); } (void)0
/*
*  Does the index interval i_:i_+n_-1 spans more than one process rows
*  or columns ?
*
*  Mspan returns 0 (false) when the data is replicated (srcproc_ < 0) or
*  when there is only one process row or column in the process grid.
*/
#define    Mspan( n_, i_, inb_, nb_, srcproc_, nprocs_ ) \
              ( ( (srcproc_) >= 0 ) && ( ( (nprocs_) > 1 ) && \
              ( ( (i_) < (inb_) ) ? \
                ( (i_) + (n_) > (inb_) ) : \
                ( (i_) + (n_) > (inb_) + \
                  ( ( (i_) - (inb_) ) / (nb_) + 1 ) * (nb_) ) ) ) )
/*
*  Mindxl2g computes the global index ig_ corresponding to the local
*  index il_ in process proc_.
*/
#define    Mindxl2g( ig_, il_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 if( (proc_) == (srcproc_) ) \
                 { \
                    if( (il_) < (inb_) ) ig_ = (il_); \
                    else                 ig_ = (il_) + \
                       (nb_)*((nprocs_)-1)*( ((il_)-(inb_)) / (nb_) + 1 ); \
                 } \
                 else if( (proc_) < (srcproc_) ) \
                 { \
                    ig_ = (il_) + (inb_) + \
                          (nb_)*(  ((nprocs_)-1)*((il_)/(nb_)) + \
                                   (proc_)-(srcproc_)-1+(nprocs_) ); \
                 } \
                 else \
                 { \
                    ig_ =  (il_) + (inb_) + \
                           (nb_)*( ((nprocs_)-1)*((il_)/(nb_)) + \
                           (proc_)-(srcproc_)-1 ); \
                 } \
              } \
              else \
              { \
                 ig_ = (il_); \
              } \
           } (void)0
/*
*  Mindxg2p returns the process coodinate owning the entry globally
*  indexed by ig_.
*/
#define    Mindxg2p( ig_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (ig_) >= (inb_) ) && ( (srcproc_) >= 0 ) && \
                  ( (nprocs_) > 1 ) ) \
              { \
                 proc_  = (srcproc_) + 1 + ( (ig_)-(inb_) ) / (nb_); \
                 proc_ -= ( proc_ / (nprocs_) ) * (nprocs_); \
              } \
              else \
              { \
                 proc_ = (srcproc_); \
              } \
           } (void)0
/*
*  Mnumroc computes the # of local indexes np_ residing in the process
*  of coordinate proc_ corresponding to the interval of global indexes
*  i_:i_+n_-1 assuming that the global index 0 resides in  the process
*  srcproc_, and that the indexes are distributed from  srcproc_ using
*  the parameters inb_, nb_ and nprocs_.
*/
#define    Mnumroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    inb__ += nblk__*(nb_); \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) \
                    { if( (proc_) == src__ ) np_ = (n_); else np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__    = n__ / (nb_) + 1; \
                       mydist__ -= nblk__ - \
                          ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < 0 ) \
                       { \
                          if( (proc_) != src__ ) \
                             np_ = (nb_) + (nb_) * quot__; \
                          else \
                             np_ = inb__ + (nb_) * quot__; \
                       } \
                       else if( mydist__ > 0 ) \
                       { \
                          np_ = (nb_) * quot__; \
                       } \
                       else \
                       { \
                          if( (proc_) != src__ ) \
                             np_ = n__ + (nb_) + (nb_) * ( quot__ - nblk__ ); \
                          else \
                             np_ = (n_) +        (nb_) * ( quot__ - nblk__ ); \
                       } \
                    } \
                 } \
                 else \
                 { \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) \
                    { if( (proc_) == (srcproc_) ) np_ = (n_); else np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__    = n__ / (nb_) + 1; \
                       mydist__ -= nblk__ - \
                          ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < 0 ) \
                       { \
                          if( (proc_) != (srcproc_) ) \
                             np_ = (nb_) + (nb_) * quot__; \
                          else \
                             np_ = inb__ + (nb_) * quot__; \
                       } \
                       else if( mydist__ > 0 ) \
                       { \
                          np_ = (nb_) * quot__; \
                       } \
                       else \
                       { \
                          if( (proc_) != (srcproc_) ) \
                             np_ = n__ + (nb_) + (nb_) * ( quot__ - nblk__ ); \
                          else \
                             np_ = (n_) +        (nb_) * ( quot__ - nblk__ ); \
                       } \
                    } \
                 } \
              } \
              else \
              { \
                 np_ = (n_); \
              } \
           } (void)0

#define    Mnpreroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, rem__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    if( (proc_) != src__ ) \
                    { \
                       inb__ += nblk__*(nb_); \
                       if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = (n_); } \
                       else \
                       { \
                          if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                             mydist__ += (nprocs_); \
                          nblk__ = n__ / (nb_) + 1; \
                          rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                          if( mydist__ <= rem__ ) \
                          { \
                             np_ = inb__ - (nb_) + \
                                   ( quot__ + 1 ) * mydist__ * (nb_); \
                          } \
                          else \
                          { \
                             np_ = (n_) + \
                                   ( mydist__ - (nprocs_) ) * quot__ * (nb_); \
                          } \
                       } \
                    } \
                    else \
                    { \
                       np_ = 0; \
                    } \
                 } \
                 else \
                 { \
                    if( (proc_) != (srcproc_) ) \
                    { \
                       if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = (n_); } \
                       else \
                       { \
                          if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                             mydist__ += (nprocs_); \
                          nblk__ = n__ / (nb_) + 1; \
                          rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                          if( mydist__ <= rem__ ) \
                          { \
                             np_ = inb__ - (nb_) + \
                                   ( quot__ + 1 ) * mydist__ * (nb_); \
                          } \
                          else \
                          { \
                             np_ = (n_) + \
                                   ( mydist__ - (nprocs_) ) * quot__ * (nb_); \
                          } \
                       } \
                    } \
                    else \
                    { \
                       np_ = 0; \
                    } \
                 } \
              } \
              else \
              { \
                 np_ = 0; \
              } \
           } (void)0

#define    Mnnxtroc( np_, n_, i_, inb_, nb_, proc_, srcproc_, nprocs_ ) \
           { \
              if( ( (srcproc_) >= 0 ) && ( (nprocs_) > 1 ) ) \
              { \
                 Int inb__, mydist__, n__, nblk__, quot__, rem__, src__; \
                 if( ( inb__ = (inb_) - (i_) ) <= 0 ) \
                 { \
                    src__  = (srcproc_) + ( nblk__ = (-inb__) / (nb_) + 1 ); \
                    src__ -= ( src__ / (nprocs_) ) * (nprocs_); \
                    inb__ += nblk__*(nb_); \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - src__ ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__ = n__ / (nb_) + 1; \
                       rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < rem__ ) \
                       { \
                          np_ = n__ - ( quot__ * mydist__ + \
                                        quot__ + mydist__ ) * (nb_); \
                       } \
                       else \
                       { \
                          np_ = ( (nprocs_) - 1 - mydist__ ) * quot__ * (nb_); \
                       } \
                    } \
                 } \
                 else \
                 { \
                    if( ( n__ = (n_) - inb__ ) <= 0 ) { np_ = 0; } \
                    else \
                    { \
                       if( ( mydist__ = (proc_) - (srcproc_) ) < 0 ) \
                          mydist__ += (nprocs_); \
                       nblk__ = n__ / (nb_) + 1; \
                       rem__  = nblk__ - \
                             ( quot__ = ( nblk__ / (nprocs_) ) ) * (nprocs_); \
                       if( mydist__ < rem__ ) \
                       { \
                          np_ = n__ - ( quot__ * mydist__ + \
                                        quot__ + mydist__ ) * (nb_); \
                       } \
                       else \
                       { \
                          np_ = ( (nprocs_) - 1 - mydist__ ) * quot__ * (nb_); \
                       } \
                    } \
                 } \
              } \
              else \
              { np_ = 0; } \
           } (void)0


#define    Minfog2l( i_, j_, desc_, nr_, nc_, r_, c_, ii_, jj_, pr_, pc_ ) \
           { \
              Int quot__, i__, imb__, inb__, j__, mb__, mydist__, \
                  nb__, nblk__, src__; \
              imb__ = desc_[IMB_]; mb__ = desc_[MB_]; pr_ = desc_[RSRC_]; \
              if( ( pr_ >= 0 ) && ( nr_ > 1 ) ) \
              { \
                 if( ( i__ = (i_) - imb__ ) < 0 ) \
                 { ii_ = ( r_ == pr_ ? (i_) : 0 ); } \
                 else \
                 { \
                    src__     = pr_; \
                    pr_      += ( nblk__ = i__ / mb__ + 1 ); \
                    pr_      -= ( pr_ / nr_ ) * nr_; \
                    if( ( mydist__ = r_ - src__ ) < 0 ) mydist__ += nr_; \
                    if( mydist__ >= nblk__ - ( quot__ = nblk__ / nr_ ) * nr_ ) \
                    { \
                       if( r_ != src__ ) ii_ =  mb__; \
                       else              ii_ = imb__; \
                       if( r_ != pr_ ) \
                          ii_ += ( quot__ - 1 ) * mb__; \
                       else \
                          ii_ += i__ + ( quot__ - nblk__ ) * mb__; \
                    } \
                    else \
                    { \
                       if( r_ != src__ ) ii_ =  mb__ + quot__ * mb__; \
                       else              ii_ = imb__ + quot__ * mb__; \
                    } \
                 } \
              } \
              else \
              { \
                 ii_ = (i_); \
              } \
              inb__ = desc_[INB_]; nb__ = desc_[NB_]; pc_ = desc_[CSRC_]; \
              if( ( pc_ >= 0 ) && ( nc_ > 1 ) ) \
              { \
                 if( ( j__ = (j_) - inb__ ) < 0 ) \
                 { jj_ = ( c_ == pc_ ? (j_) : 0 ); } \
                 else \
                 { \
                    src__     = pc_; \
                    pc_      += ( nblk__ = j__ / nb__ + 1 ); \
                    pc_      -= ( pc_ / nc_ ) * nc_; \
                    if( ( mydist__ = c_ - src__ ) < 0 ) mydist__ += nc_; \
                    if( mydist__ >= nblk__ - ( quot__ = nblk__ / nc_ ) * nc_ ) \
                    { \
                       if( c_ != src__ ) jj_ =  nb__; \
                       else              jj_ = inb__; \
                       if( c_ != pc_ ) \
                          jj_ += ( quot__ - 1 ) * nb__; \
                       else \
                          jj_ += j__ + ( quot__ - nblk__ ) * nb__; \
                    } \
                    else \
                    { \
                       if( c_ != src__ ) jj_ =  nb__ + quot__ * nb__; \
                       else              jj_ = inb__ + quot__ * nb__; \
                    } \
                 } \
              } \
              else \
              { \
                 jj_ = (j_); \
              } \
           } (void)0

/*
*  The following macros initialize or translate descriptors.
*/
#define    MDescSet( desc, m, n, imb, inb, mb, nb, rsrc, csrc, ictxt, lld ) \
           { \
              (desc)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
              (desc)[CTXT_ ] = (ictxt); \
              (desc)[M_    ] = (m);     \
              (desc)[N_    ] = (n);     \
              (desc)[IMB_  ] = (imb);   \
              (desc)[INB_  ] = (inb);   \
              (desc)[MB_   ] = (mb);    \
              (desc)[NB_   ] = (nb);    \
              (desc)[RSRC_ ] = (rsrc);  \
              (desc)[CSRC_ ] = (csrc);  \
              (desc)[LLD_  ] = (lld);   \
           } (void)0

#define    MDescCopy(DescIn, DescOut) \
           { \
              (DescOut)[DTYPE_] = (DescIn)[DTYPE_];    \
              (DescOut)[M_    ] = (DescIn)[M_    ];    \
              (DescOut)[N_    ] = (DescIn)[N_    ];    \
              (DescOut)[IMB_  ] = (DescIn)[IMB_  ];    \
              (DescOut)[INB_  ] = (DescIn)[INB_  ];    \
              (DescOut)[MB_   ] = (DescIn)[MB_   ];    \
              (DescOut)[NB_   ] = (DescIn)[NB_   ];    \
              (DescOut)[RSRC_ ] = (DescIn)[RSRC_ ];    \
              (DescOut)[CSRC_ ] = (DescIn)[CSRC_ ];    \
              (DescOut)[CTXT_ ] = (DescIn)[CTXT_ ];    \
              (DescOut)[LLD_  ] = (DescIn)[LLD_  ];    \
           } (void)0

#define    MDescTrans(DescIn, DescOut) \
           { \
              if ( (DescIn)[DTYPE_] == BLOCK_CYCLIC_2D ) \
              { \
                 (DescOut)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
                 (DescOut)[M_    ] = (DescIn)[M1_    ];   \
                 (DescOut)[N_    ] = (DescIn)[N1_    ];   \
                 (DescOut)[IMB_  ] = (DescIn)[MB1_   ];   \
                 (DescOut)[INB_  ] = (DescIn)[NB1_   ];   \
                 (DescOut)[MB_   ] = (DescIn)[MB1_   ];   \
                 (DescOut)[NB_   ] = (DescIn)[NB1_   ];   \
                 (DescOut)[RSRC_ ] = (DescIn)[RSRC1_ ];   \
                 (DescOut)[CSRC_ ] = (DescIn)[CSRC1_ ];   \
                 (DescOut)[CTXT_ ] = (DescIn)[CTXT1_ ];   \
                 (DescOut)[LLD_  ] = (DescIn)[LLD1_  ];   \
              } \
              else if ( (DescIn)[DTYPE_] == BLOCK_CYCLIC_2D_INB ) \
              { \
                 (DescOut)[DTYPE_] = BLOCK_CYCLIC_2D_INB; \
                 (DescOut)[M_    ] = (DescIn)[M_    ];    \
                 (DescOut)[N_    ] = (DescIn)[N_    ];    \
                 (DescOut)[IMB_  ] = (DescIn)[IMB_  ];    \
                 (DescOut)[INB_  ] = (DescIn)[INB_  ];    \
                 (DescOut)[MB_   ] = (DescIn)[MB_   ];    \
                 (DescOut)[NB_   ] = (DescIn)[NB_   ];    \
                 (DescOut)[RSRC_ ] = (DescIn)[RSRC_ ];    \
                 (DescOut)[CSRC_ ] = (DescIn)[CSRC_ ];    \
                 (DescOut)[CTXT_ ] = (DescIn)[CTXT_ ];    \
                 (DescOut)[LLD_  ] = (DescIn)[LLD_  ];    \
              } \
              else \
              { \
                 (DescOut)[DTYPE_] = (DescIn)[0]; \
                 (DescOut)[CTXT_ ] = (DescIn)[1]; \
                 (DescOut)[M_    ] = 0;           \
                 (DescOut)[N_    ] = 0;           \
                 (DescOut)[IMB_  ] = 1;           \
                 (DescOut)[INB_  ] = 1;           \
                 (DescOut)[MB_   ] = 1;           \
                 (DescOut)[NB_   ] = 1;           \
                 (DescOut)[RSRC_ ] = 0;           \
                 (DescOut)[CSRC_ ] = 0;           \
                 (DescOut)[LLD_  ] = 1;           \
              } \
           } (void)0

#define    MIndxTrans( I, J, i, j ) \
           { \
              i = *I - 1; \
              j = *J - 1; \
           } (void)0

#if( _F2C_CALL_ == _F2C_ADD_ )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine. No redefinition is necessary  to  have
*  the following FORTRAN to C interface:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo_(...)
*
*  This is the PBLAS default.
*/

#endif

#if( _F2C_CALL_ == _F2C_F77ISF2C )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine for systems where  the FORTRAN compiler
*  is actually f2c (a FORTRAN to C conversion utility).
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo__(...)
*/

#endif

#if( _F2C_CALL_ == _F2C_UPCASE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          PDFOO(...)
*/
#define    immadd_             IMMADD
#define    smmadd_             SMMADD
#define    emmadd_             EMMADD
#define    dmmadd_             DMMADD
#define    emmadd_             EMMADD
#define    cmmadd_             CMMADD
#define    ymmadd_             YMMADD
#define    zmmadd_             ZMMADD
#define    ymmadd_             YMMADD

#define    immtadd_            IMMTADD
#define    smmtadd_            SMMTADD
#define    emmtadd_            EMMTADD
#define    dmmtadd_            DMMTADD
#define    emmtadd_            EMMTADD
#define    cmmtadd_            CMMTADD
#define    ymmtadd_            YMMTADD
#define    zmmtadd_            ZMMTADD
#define    ymmtadd_            YMMTADD

#define    smmcadd_            SMMCADD
#define    emmcadd_            EMMCADD
#define    dmmcadd_            DMMCADD
#define    emmcadd_            EMMCADD
#define    cmmcadd_            CMMCADD
#define    ymmcadd_            YMMCADD
#define    zmmcadd_            ZMMCADD
#define    ymmcadd_            YMMCADD

#define    smmtcadd_           SMMTCADD
#define    emmtcadd_           EMMTCADD
#define    dmmtcadd_           DMMTCADD
#define    emmtcadd_           EMMTCADD
#define    cmmtcadd_           CMMTCADD
#define    ymmtcadd_           YMMTCADD
#define    zmmtcadd_           ZMMTCADD
#define    ymmtcadd_           YMMTCADD

#define    immdda_             IMMDDA
#define    smmdda_             SMMDDA
#define    emmdda_             EMMDDA
#define    dmmdda_             DMMDDA
#define    emmdda_             EMMDDA
#define    cmmdda_             CMMDDA
#define    ymmdda_             YMMDDA
#define    zmmdda_             ZMMDDA
#define    ymmdda_             YMMDDA

#define    smmddac_            SMMDDAC
#define    emmddac_            EMMDDAC
#define    dmmddac_            DMMDDAC
#define    emmddac_            EMMDDAC
#define    cmmddac_            CMMDDAC
#define    ymmddac_            YMMDDAC
#define    zmmddac_            ZMMDDAC
#define    ymmddac_            YMMDDAC

#define    immddat_            IMMDDAT
#define    smmddat_            SMMDDAT
#define    emmddat_            EMMDDAT
#define    dmmddat_            DMMDDAT
#define    emmddat_            EMMDDAT
#define    cmmddat_            CMMDDAT
#define    ymmddat_            YMMDDAT
#define    zmmddat_            ZMMDDAT
#define    ymmddat_            YMMDDAT

#define    smmddact_           SMMDDACT
#define    emmddact_           EMMDDACT
#define    dmmddact_           DMMDDACT
#define    emmddact_           EMMDDACT
#define    cmmddact_           CMMDDACT
#define    ymmddact_           YMMDDACT
#define    zmmddact_           ZMMDDACT
#define    ymmddact_           YMMDDACT

#define    sasqrtb_            SASQRTB
#define    easqrtb_            EASQRTB
#define    dasqrtb_            DASQRTB
#define    easqrtb_            EASQRTB

#define    sset_               SSET
#define    eset_               ESET
#define    dset_               DSET
#define    eset_               ESET
#define    cset_               CSET
#define    yset_               YSET
#define    zset_               ZSET
#define    yset_               YSET

#define    svasum_             SVASUM
#define    evasum_             EVASUM
#define    dvasum_             DVASUM
#define    evasum_             EVASUM
#define    scvasum_            SCVASUM
#define    eyvasum_            EYVASUM
#define    dzvasum_            DZVASUM
#define    eyvasum_            EYVASUM

#define    sascal_             SASCAL
#define    eascal_             EASCAL
#define    dascal_             DASCAL
#define    eascal_             EASCAL

#define    scshft_             SCSHFT
#define    ecshft_             ECSHFT
#define    dcshft_             DCSHFT
#define    ecshft_             ECSHFT
#define    ccshft_             CCSHFT
#define    ycshft_             YCSHFT
#define    zcshft_             ZCSHFT
#define    ycshft_             YCSHFT

#define    srshft_             SRSHFT
#define    ershft_             ERSHFT
#define    drshft_             DRSHFT
#define    ershft_             ERSHFT
#define    crshft_             CRSHFT
#define    yrshft_             YRSHFT
#define    zrshft_             ZRSHFT
#define    yrshft_             YRSHFT

#define    svvdot_             SVVDOT
#define    evvdot_             EVVDOT
#define    dvvdot_             DVVDOT
#define    evvdot_             EVVDOT
#define    cvvdotc_            CVVDOTC
#define    yvvdotc_            YVVDOTC
#define    cvvdotu_            CVVDOTU
#define    yvvdotu_            YVVDOTU
#define    zvvdotc_            ZVVDOTC
#define    yvvdotc_            YVVDOTC
#define    zvvdotu_            ZVVDOTU
#define    yvvdotu_            YVVDOTU

#define    stzpad_             STZPAD
#define    etzpad_             ETZPAD
#define    dtzpad_             DTZPAD
#define    etzpad_             ETZPAD
#define    ctzpad_             CTZPAD
#define    ytzpad_             YTZPAD
#define    ztzpad_             ZTZPAD
#define    ytzpad_             YTZPAD

#define    stzpadcpy_          STZPADCPY
#define    etzpadcpy_          ETZPADCPY
#define    dtzpadcpy_          DTZPADCPY
#define    etzpadcpy_          ETZPADCPY
#define    ctzpadcpy_          CTZPADCPY
#define    ytzpadcpy_          YTZPADCPY
#define    ztzpadcpy_          ZTZPADCPY
#define    ytzpadcpy_          YTZPADCPY

#define    stzscal_            STZSCAL
#define    etzscal_            ETZSCAL
#define    dtzscal_            DTZSCAL
#define    etzscal_            ETZSCAL
#define    ctzscal_            CTZSCAL
#define    ytzscal_            YTZSCAL
#define    ztzscal_            ZTZSCAL
#define    ytzscal_            YTZSCAL

#define    chescal_            CHESCAL
#define    yhescal_            YHESCAL
#define    zhescal_            ZHESCAL
#define    yhescal_            YHESCAL

#define    ctzcnjg_            CTZCNJG
#define    ytzcnjg_            YTZCNJG
#define    ztzcnjg_            ZTZCNJG
#define    ytzcnjg_            YTZCNJG

#define    sagemv_             SAGEMV
#define    eagemv_             EAGEMV
#define    dagemv_             DAGEMV
#define    eagemv_             EAGEMV
#define    cagemv_             CAGEMV
#define    yagemv_             YAGEMV
#define    zagemv_             ZAGEMV
#define    yagemv_             YAGEMV

#define    sasymv_             SASYMV
#define    easymv_             EASYMV
#define    dasymv_             DASYMV
#define    easymv_             EASYMV
#define    casymv_             CASYMV
#define    yasymv_             YASYMV
#define    zasymv_             ZASYMV
#define    yasymv_             YASYMV
#define    cahemv_             CAHEMV
#define    yahemv_             YAHEMV
#define    zahemv_             ZAHEMV
#define    yahemv_             YAHEMV

#define    satrmv_             SATRMV
#define    eatrmv_             EATRMV
#define    datrmv_             DATRMV
#define    eatrmv_             EATRMV
#define    catrmv_             CATRMV
#define    yatrmv_             YATRMV
#define    zatrmv_             ZATRMV
#define    yatrmv_             YATRMV

#define    csymv_              CSYMV
#define    ysymv_              YSYMV
#define    zsymv_              ZSYMV
#define    ysymv_              YSYMV

#define    csyr_               CSYR
#define    ysyr_               YSYR
#define    zsyr_               ZSYR
#define    ysyr_               YSYR

#define    csyr2_              CSYR2
#define    ysyr2_              YSYR2
#define    zsyr2_              ZSYR2
#define    ysyr2_              YSYR2

#endif

#if( _F2C_CALL_ == _F2C_NOCHANGE )
/*
*  These defines  set  up  the  naming scheme required to have a FORTRAN
*  routine called by a C routine with the following  FORTRAN to C inter-
*  face:
*
*           FORTRAN DECLARATION            C CALL
*           SUBROUTINE PDFOO(...)          pdfoo(...)
*/
#define    immadd_             immadd
#define    smmadd_             smmadd
#define    emmadd_             emmadd
#define    dmmadd_             dmmadd
#define    emmadd_             emmadd
#define    cmmadd_             cmmadd
#define    ymmadd_             ymmadd
#define    zmmadd_             zmmadd
#define    ymmadd_             ymmadd

#define    immtadd_            immtadd
#define    smmtadd_            smmtadd
#define    emmtadd_            emmtadd
#define    dmmtadd_            dmmtadd
#define    emmtadd_            emmtadd
#define    cmmtadd_            cmmtadd
#define    ymmtadd_            ymmtadd
#define    zmmtadd_            zmmtadd
#define    ymmtadd_            ymmtadd

#define    smmcadd_            smmcadd
#define    emmcadd_            emmcadd
#define    dmmcadd_            dmmcadd
#define    emmcadd_            emmcadd
#define    cmmcadd_            cmmcadd
#define    ymmcadd_            ymmcadd
#define    zmmcadd_            zmmcadd
#define    ymmcadd_            ymmcadd

#define    smmtcadd_           smmtcadd
#define    emmtcadd_           emmtcadd
#define    dmmtcadd_           dmmtcadd
#define    emmtcadd_           emmtcadd
#define    cmmtcadd_           cmmtcadd
#define    ymmtcadd_           ymmtcadd
#define    zmmtcadd_           zmmtcadd
#define    ymmtcadd_           ymmtcadd

#define    immdda_             immdda
#define    smmdda_             smmdda
#define    emmdda_             emmdda
#define    dmmdda_             dmmdda
#define    emmdda_             emmdda
#define    cmmdda_             cmmdda
#define    ymmdda_             ymmdda
#define    zmmdda_             zmmdda
#define    ymmdda_             ymmdda

#define    smmddac_            smmddac
#define    emmddac_            emmddac
#define    dmmddac_            dmmddac
#define    emmddac_            emmddac
#define    cmmddac_            cmmddac
#define    ymmddac_            ymmddac
#define    zmmddac_            zmmddac
#define    ymmddac_            ymmddac

#define    immddat_            immddat
#define    smmddat_            smmddat
#define    emmddat_            emmddat
#define    dmmddat_            dmmddat
#define    emmddat_            emmddat
#define    cmmddat_            cmmddat
#define    ymmddat_            ymmddat
#define    zmmddat_            zmmddat
#define    ymmddat_            ymmddat

#define    smmddact_           smmddact
#define    emmddact_           emmddact
#define    dmmddact_           dmmddact
#define    emmddact_           emmddact
#define    cmmddact_           cmmddact
#define    ymmddact_           ymmddact
#define    zmmddact_           zmmddact
#define    ymmddact_           ymmddact

#define    sasqrtb_            sasqrtb
#define    easqrtb_            easqrtb
#define    dasqrtb_            dasqrtb
#define    easqrtb_            easqrtb

#define    sset_               sset
#define    eset_               eset
#define    dset_               dset
#define    eset_               eset
#define    cset_               cset
#define    yset_               yset
#define    zset_               zset
#define    yset_               yset

#define    svasum_             svasum
#define    evasum_             evasum
#define    dvasum_             dvasum
#define    evasum_             evasum
#define    scvasum_            scvasum
#define    eyvasum_            eyvasum
#define    dzvasum_            dzvasum
#define    eyvasum_            eyvasum

#define    sascal_             sascal
#define    eascal_             eascal
#define    dascal_             dascal
#define    eascal_             eascal

#define    scshft_             scshft
#define    ecshft_             ecshft
#define    dcshft_             dcshft
#define    ecshft_             ecshft
#define    ccshft_             ccshft
#define    ycshft_             ycshft
#define    zcshft_             zcshft
#define    ycshft_             ycshft

#define    srshft_             srshft
#define    ershft_             ershft
#define    drshft_             drshft
#define    ershft_             ershft
#define    crshft_             crshft
#define    yrshft_             yrshft
#define    zrshft_             zrshft
#define    yrshft_             yrshft

#define    svvdot_             svvdot
#define    evvdot_             evvdot
#define    dvvdot_             dvvdot
#define    evvdot_             evvdot
#define    cvvdotc_            cvvdotc
#define    yvvdotc_            yvvdotc
#define    cvvdotu_            cvvdotu
#define    yvvdotu_            yvvdotu
#define    zvvdotc_            zvvdotc
#define    yvvdotc_            yvvdotc
#define    zvvdotu_            zvvdotu
#define    yvvdotu_            yvvdotu

#define    stzpad_             stzpad
#define    etzpad_             etzpad
#define    dtzpad_             dtzpad
#define    etzpad_             etzpad
#define    ctzpad_             ctzpad
#define    ytzpad_             ytzpad
#define    ztzpad_             ztzpad
#define    ytzpad_             ytzpad

#define    stzpadcpy_          stzpadcpy
#define    etzpadcpy_          etzpadcpy
#define    dtzpadcpy_          dtzpadcpy
#define    etzpadcpy_          etzpadcpy
#define    ctzpadcpy_          ctzpadcpy
#define    ytzpadcpy_          ytzpadcpy
#define    ztzpadcpy_          ztzpadcpy
#define    ytzpadcpy_          ytzpadcpy

#define    stzscal_            stzscal
#define    etzscal_            etzscal
#define    dtzscal_            dtzscal
#define    etzscal_            etzscal
#define    ctzscal_            ctzscal
#define    ytzscal_            ytzscal
#define    ztzscal_            ztzscal
#define    ytzscal_            ytzscal

#define    chescal_            chescal
#define    yhescal_            yhescal
#define    zhescal_            zhescal
#define    yhescal_            yhescal

#define    ctzcnjg_            ctzcnjg
#define    ytzcnjg_            ytzcnjg
#define    ztzcnjg_            ztzcnjg
#define    ytzcnjg_            ytzcnjg

#define    sagemv_             sagemv
#define    eagemv_             eagemv
#define    dagemv_             dagemv
#define    eagemv_             eagemv
#define    cagemv_             cagemv
#define    yagemv_             yagemv
#define    zagemv_             zagemv
#define    yagemv_             yagemv

#define    sasymv_             sasymv
#define    easymv_             easymv
#define    dasymv_             dasymv
#define    easymv_             easymv
#define    casymv_             casymv
#define    yasymv_             yasymv
#define    zasymv_             zasymv
#define    yasymv_             yasymv
#define    cahemv_             cahemv
#define    yahemv_             yahemv
#define    zahemv_             zahemv
#define    yahemv_             yahemv

#define    satrmv_             satrmv
#define    eatrmv_             eatrmv
#define    datrmv_             datrmv
#define    eatrmv_             eatrmv
#define    catrmv_             catrmv
#define    yatrmv_             yatrmv
#define    zatrmv_             zatrmv
#define    yatrmv_             yatrmv

#define    csymv_              csymv
#define    ysymv_              ysymv
#define    zsymv_              zsymv
#define    ysymv_              ysymv

#define    csyr_               csyr
#define    ysyr_               ysyr
#define    zsyr_               zsyr
#define    ysyr_               ysyr

#define    csyr2_              csyr2
#define    ysyr2_              ysyr2
#define    zsyr2_              zsyr2
#define    ysyr2_              ysyr2

#endif
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

F_VOID_FCT     immadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmadd_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmcadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmtadd_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmtcadd_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmdda_         ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddac_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     immddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     smmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddat_        ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     smmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     emmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zmmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     ymmddact_       ( Int *,     Int *,     char *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     sasqrtb_        ( float *,   float *,   float * );
F_VOID_FCT     easqrtb_        ( EREAL *,   EREAL *,   EREAL * );
F_VOID_FCT     dasqrtb_        ( double *,  double *,  double * );
F_VOID_FCT     easqrtb_        ( EREAL *,  EREAL *,  EREAL * );

F_VOID_FCT     sset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zset_           ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yset_           ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     svasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     evasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dvasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     evasum_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     scvasum_        ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eyvasum_        ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dzvasum_        ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eyvasum_        ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sascal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eascal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dascal_         ( Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     eascal_         ( Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     scshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ecshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     dcshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ecshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ccshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ycshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     zcshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ycshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );

F_VOID_FCT     srshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ershft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     drshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     ershft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     crshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     yrshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     zrshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );
F_VOID_FCT     yrshft_         ( Int *,     Int *,     Int *,
                                 char *,    Int * );

F_VOID_FCT     svvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     evvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dvvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     evvdot_         ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     cvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yvvdotu_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     yvvdotc_        ( Int *,     char *,    char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     stzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     etzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     dtzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     etzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ctzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ytzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ztzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );
F_VOID_FCT     ytzpad_         ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 char *,    char *,    Int * );

F_VOID_FCT     stzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     etzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     dtzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     etzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ctzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ytzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ztzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ytzpadcpy_      ( F_CHAR_T,  F_CHAR_T,  Int *,
                                 Int *,     Int *,     char *,
                                 Int *,     char *,    Int * );

F_VOID_FCT     stzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     etzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dtzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     etzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ctzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ytzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ztzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ytzscal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     chescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yhescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zhescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yhescal_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     ctzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ytzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ztzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ytzcnjg_        ( F_CHAR_T,  Int *,     Int *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     sagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     eagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     dagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     eagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     cagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     yagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     zagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );
F_VOID_FCT     yagemv_         ( F_CHAR_T,  Int *,     Int *,
                                 char *,    char *,    Int *,
                                 char *,    Int *,     char *,
                                 char *,    Int * );

F_VOID_FCT     sasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     easymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     dasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     easymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     casymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yasymv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     cahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     yahemv_         ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     satrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     eatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     datrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     eatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     catrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     zatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );
F_VOID_FCT     yatrmv_         ( F_CHAR_T,  F_CHAR_T,  F_CHAR_T,
                                 Int *,     char *,    char *,
                                 Int *,     char *,    Int *,
                                 char *,    char *,    Int * );

F_VOID_FCT     csymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ysymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     zsymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );
F_VOID_FCT     ysymv_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    char *,
                                 Int * );

F_VOID_FCT     csyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     ysyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     zsyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );
F_VOID_FCT     ysyr_           ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int * );

F_VOID_FCT     csyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ysyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     zsyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );
F_VOID_FCT     ysyr2_          ( F_CHAR_T,  Int *,     char *,
                                 char *,    Int *,     char *,
                                 Int *,     char *,    Int * );

void           PB_Ctzsyr       ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzher       ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsyr2      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctzher2      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctztrmv      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzatrmv     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsymv      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzhemv      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzasymv     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzahemv     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );

void           PB_Ctzsyrk      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzherk      ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsyr2k     ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctzher2k     ( PBTYP_T *, char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 char *,    Int );
void           PB_Ctztrmm      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       Int,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
void           PB_Ctzsymm      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );
void           PB_Ctzhemm      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int );

void           PB_CpswapNN     ( PBTYP_T *, Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     Int );
void           PB_CpswapND     ( PBTYP_T *, Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     Int );
void           PB_Cpdot11      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpdotNN      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpdotND      ( PBTYP_T *, Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 Int,       VVDOT_T );
void           PB_CpaxpbyNN    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_CpaxpbyND    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_CpaxpbyDN    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );
void           PB_Cpaxpby      ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    char *,
                                 Int,       Int,       Int *,
                                 char * );

void           PB_Cpsyr        ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       Int,
                                 Int *,     TZSYR_T );
void           PB_Cpsyr2       ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     TZSYR2_T );
void           PB_Cptrm        ( PBTYP_T *, PBTYP_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int,       TZTRM_T );
void           PB_Cpsym        ( PBTYP_T *, PBTYP_T *, char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 Int,       char *,    Int,
                                 char *,    Int,       char *,
                                 Int,       TZSYM_T );
void           PB_Cpgeadd      ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cptradd      ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cptran       ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cptrsv       ( PBTYP_T *, Int,       char *,
                                 char *,    char *,    Int,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int );
void           PB_Cptrsm       ( PBTYP_T *, Int,       char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 char *,    Int );

void           PB_CpgemmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpgemmAC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpgemmBC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_CpsymmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CpsymmBC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CpsyrkA      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_CpsyrkAC     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cpsyr2kA     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_Cpsyr2kAC    ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    char *,
                                 Int,       Int,       Int * );
void           PB_CptrmmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrmmB      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrsmAB     ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
void           PB_CptrsmAB0    ( PBTYP_T *, char *,    char *,
                                 char *,    Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 Int,       Int,       Int *,
                                 char * *,  Int *,     Int * );
void           PB_CptrsmAB1    ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    Int,
                                 Int,       char *,    char *,
                                 Int,       Int,       Int *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int * );
void           PB_CptrsmB      ( PBTYP_T *, char *,    char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int *,     char *,    Int,
                                 Int,       Int * );
#else

F_VOID_FCT     immadd_         ();
F_VOID_FCT     smmadd_         ();
F_VOID_FCT     emmadd_         ();
F_VOID_FCT     dmmadd_         ();
F_VOID_FCT     emmadd_         ();
F_VOID_FCT     cmmadd_         ();
F_VOID_FCT     ymmadd_         ();
F_VOID_FCT     zmmadd_         ();
F_VOID_FCT     ymmadd_         ();

F_VOID_FCT     smmcadd_        ();
F_VOID_FCT     emmcadd_        ();
F_VOID_FCT     dmmcadd_        ();
F_VOID_FCT     emmcadd_        ();
F_VOID_FCT     cmmcadd_        ();
F_VOID_FCT     ymmcadd_        ();
F_VOID_FCT     zmmcadd_        ();
F_VOID_FCT     ymmcadd_        ();

F_VOID_FCT     immtadd_        ();
F_VOID_FCT     smmtadd_        ();
F_VOID_FCT     emmtadd_        ();
F_VOID_FCT     dmmtadd_        ();
F_VOID_FCT     emmtadd_        ();
F_VOID_FCT     cmmtadd_        ();
F_VOID_FCT     ymmtadd_        ();
F_VOID_FCT     zmmtadd_        ();
F_VOID_FCT     ymmtadd_        ();

F_VOID_FCT     smmtcadd_       ();
F_VOID_FCT     emmtcadd_       ();
F_VOID_FCT     dmmtcadd_       ();
F_VOID_FCT     emmtcadd_       ();
F_VOID_FCT     cmmtcadd_       ();
F_VOID_FCT     ymmtcadd_       ();
F_VOID_FCT     zmmtcadd_       ();
F_VOID_FCT     ymmtcadd_       ();

F_VOID_FCT     immdda_         ();
F_VOID_FCT     smmdda_         ();
F_VOID_FCT     emmdda_         ();
F_VOID_FCT     dmmdda_         ();
F_VOID_FCT     emmdda_         ();
F_VOID_FCT     cmmdda_         ();
F_VOID_FCT     ymmdda_         ();
F_VOID_FCT     zmmdda_         ();
F_VOID_FCT     ymmdda_         ();

F_VOID_FCT     smmddac_        ();
F_VOID_FCT     emmddac_        ();
F_VOID_FCT     dmmddac_        ();
F_VOID_FCT     emmddac_        ();
F_VOID_FCT     cmmddac_        ();
F_VOID_FCT     ymmddac_        ();
F_VOID_FCT     zmmddac_        ();
F_VOID_FCT     ymmddac_        ();

F_VOID_FCT     immddat_        ();
F_VOID_FCT     smmddat_        ();
F_VOID_FCT     emmddat_        ();
F_VOID_FCT     dmmddat_        ();
F_VOID_FCT     emmddat_        ();
F_VOID_FCT     cmmddat_        ();
F_VOID_FCT     ymmddat_        ();
F_VOID_FCT     zmmddat_        ();
F_VOID_FCT     ymmddat_        ();

F_VOID_FCT     smmddact_       ();
F_VOID_FCT     emmddact_       ();
F_VOID_FCT     dmmddact_       ();
F_VOID_FCT     emmddact_       ();
F_VOID_FCT     cmmddact_       ();
F_VOID_FCT     ymmddact_       ();
F_VOID_FCT     zmmddact_       ();
F_VOID_FCT     ymmddact_       ();

F_VOID_FCT     sasqrtb_        ();
F_VOID_FCT     easqrtb_        ();
F_VOID_FCT     dasqrtb_        ();
F_VOID_FCT     easqrtb_        ();

F_VOID_FCT     sset_           ();
F_VOID_FCT     eset_           ();
F_VOID_FCT     dset_           ();
F_VOID_FCT     eset_           ();
F_VOID_FCT     cset_           ();
F_VOID_FCT     yset_           ();
F_VOID_FCT     zset_           ();
F_VOID_FCT     yset_           ();

F_VOID_FCT     svasum_         ();
F_VOID_FCT     evasum_         ();
F_VOID_FCT     dvasum_         ();
F_VOID_FCT     evasum_         ();
F_VOID_FCT     scvasum_        ();
F_VOID_FCT     eyvasum_        ();
F_VOID_FCT     dzvasum_        ();
F_VOID_FCT     eyvasum_        ();

F_VOID_FCT     sascal_         ();
F_VOID_FCT     eascal_         ();
F_VOID_FCT     dascal_         ();
F_VOID_FCT     eascal_         ();

F_VOID_FCT     scshft_         ();
F_VOID_FCT     ecshft_         ();
F_VOID_FCT     dcshft_         ();
F_VOID_FCT     ecshft_         ();
F_VOID_FCT     ccshft_         ();
F_VOID_FCT     ycshft_         ();
F_VOID_FCT     zcshft_         ();
F_VOID_FCT     ycshft_         ();

F_VOID_FCT     srshft_         ();
F_VOID_FCT     ershft_         ();
F_VOID_FCT     drshft_         ();
F_VOID_FCT     ershft_         ();
F_VOID_FCT     crshft_         ();
F_VOID_FCT     yrshft_         ();
F_VOID_FCT     zrshft_         ();
F_VOID_FCT     yrshft_         ();

F_VOID_FCT     svvdot_         ();
F_VOID_FCT     evvdot_         ();
F_VOID_FCT     dvvdot_         ();
F_VOID_FCT     evvdot_         ();
F_VOID_FCT     cvvdotc_        ();
F_VOID_FCT     yvvdotc_        ();
F_VOID_FCT     cvvdotu_        ();
F_VOID_FCT     yvvdotu_        ();
F_VOID_FCT     zvvdotc_        ();
F_VOID_FCT     yvvdotc_        ();
F_VOID_FCT     zvvdotu_        ();
F_VOID_FCT     yvvdotu_        ();

F_VOID_FCT     stzpad_         ();
F_VOID_FCT     etzpad_         ();
F_VOID_FCT     dtzpad_         ();
F_VOID_FCT     etzpad_         ();
F_VOID_FCT     ctzpad_         ();
F_VOID_FCT     ytzpad_         ();
F_VOID_FCT     ztzpad_         ();
F_VOID_FCT     ytzpad_         ();

F_VOID_FCT     stzpadcpy_      ();
F_VOID_FCT     etzpadcpy_      ();
F_VOID_FCT     dtzpadcpy_      ();
F_VOID_FCT     etzpadcpy_      ();
F_VOID_FCT     ctzpadcpy_      ();
F_VOID_FCT     ytzpadcpy_      ();
F_VOID_FCT     ztzpadcpy_      ();
F_VOID_FCT     ytzpadcpy_      ();

F_VOID_FCT     stzscal_        ();
F_VOID_FCT     etzscal_        ();
F_VOID_FCT     dtzscal_        ();
F_VOID_FCT     etzscal_        ();
F_VOID_FCT     ctzscal_        ();
F_VOID_FCT     ytzscal_        ();
F_VOID_FCT     ztzscal_        ();
F_VOID_FCT     ytzscal_        ();

F_VOID_FCT     chescal_        ();
F_VOID_FCT     yhescal_        ();
F_VOID_FCT     zhescal_        ();
F_VOID_FCT     yhescal_        ();

F_VOID_FCT     ctzcnjg_        ();
F_VOID_FCT     ytzcnjg_        ();
F_VOID_FCT     ztzcnjg_        ();
F_VOID_FCT     ytzcnjg_        ();

F_VOID_FCT     sagemv_         ();
F_VOID_FCT     eagemv_         ();
F_VOID_FCT     dagemv_         ();
F_VOID_FCT     eagemv_         ();
F_VOID_FCT     cagemv_         ();
F_VOID_FCT     yagemv_         ();
F_VOID_FCT     zagemv_         ();
F_VOID_FCT     yagemv_         ();

F_VOID_FCT     sasymv_         ();
F_VOID_FCT     easymv_         ();
F_VOID_FCT     dasymv_         ();
F_VOID_FCT     easymv_         ();
F_VOID_FCT     casymv_         ();
F_VOID_FCT     yasymv_         ();
F_VOID_FCT     zasymv_         ();
F_VOID_FCT     yasymv_         ();
F_VOID_FCT     cahemv_         ();
F_VOID_FCT     yahemv_         ();
F_VOID_FCT     zahemv_         ();
F_VOID_FCT     yahemv_         ();

F_VOID_FCT     satrmv_         ();
F_VOID_FCT     eatrmv_         ();
F_VOID_FCT     datrmv_         ();
F_VOID_FCT     eatrmv_         ();
F_VOID_FCT     catrmv_         ();
F_VOID_FCT     yatrmv_         ();
F_VOID_FCT     zatrmv_         ();
F_VOID_FCT     yatrmv_         ();

F_VOID_FCT     csymv_          ();
F_VOID_FCT     ysymv_          ();
F_VOID_FCT     zsymv_          ();
F_VOID_FCT     ysymv_          ();

F_VOID_FCT     csyr_           ();
F_VOID_FCT     ysyr_           ();
F_VOID_FCT     zsyr_           ();
F_VOID_FCT     ysyr_           ();

F_VOID_FCT     csyr2_          ();
F_VOID_FCT     ysyr2_          ();
F_VOID_FCT     zsyr2_          ();
F_VOID_FCT     ysyr2_          ();

void           PB_Ctzsyr       ();
void           PB_Ctzher       ();
void           PB_Ctzsyr2      ();
void           PB_Ctzher2      ();
void           PB_Ctztrmv      ();
void           PB_Ctzatrmv     ();
void           PB_Ctzsymv      ();
void           PB_Ctzhemv      ();
void           PB_Ctzasymv     ();
void           PB_Ctzahemv     ();
void           PB_Ctzsyrk      ();
void           PB_Ctzherk      ();
void           PB_Ctzsyr2k     ();
void           PB_Ctzher2k     ();
void           PB_Ctztrmm      ();
void           PB_Ctzsymm      ();
void           PB_Ctzhemm      ();

void           PB_CpswapNN     ();
void           PB_CpswapND     ();
void           PB_Cpdot11      ();
void           PB_CpdotNN      ();
void           PB_CpdotND      ();
void           PB_CpaxpbyNN    ();
void           PB_CpaxpbyND    ();
void           PB_CpaxpbyDN    ();
void           PB_Cpaxpby      ();

void           PB_Cpsyr        ();
void           PB_Cpsyr2       ();
void           PB_Cptrm        ();
void           PB_Cpsym        ();
void           PB_Cpgeadd      ();
void           PB_Cptradd      ();
void           PB_Cptran       ();
void           PB_Cptrsv       ();
void           PB_Cptrsm       ();

void           PB_CpgemmAB     ();
void           PB_CpgemmAC     ();
void           PB_CpgemmBC     ();
void           PB_CpsymmAB     ();
void           PB_CpsymmBC     ();
void           PB_CpsyrkA      ();
void           PB_CpsyrkAC     ();
void           PB_Cpsyr2kA     ();
void           PB_Cpsyr2kAC    ();
void           PB_CptrmmAB     ();
void           PB_CptrmmB      ();
void           PB_CptrsmAB     ();
void           PB_CptrsmAB0    ();
void           PB_CptrsmAB1    ();
void           PB_CptrsmB      ();

#endif
                                                             /* TOOLS */
#ifdef __STDC__

Int            PB_Cgcd         ( Int,       Int );
Int            PB_Clcm         ( Int,       Int );

void           PB_Cdescset     ( Int *,     Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int );
void           PB_Cdescribe    ( Int,       Int,       Int,
                                 Int,       Int *,     Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           PB_CargFtoC     ( Int,       Int,       Int *,
                                 Int *,     Int *,     Int * );
Int            PB_Cfirstnb     ( Int,       Int,       Int,
                                 Int );
Int            PB_Clastnb      ( Int,       Int,       Int,
                                 Int );
Int            PB_Cspan        ( Int,       Int,       Int,
                                 Int,       Int,       Int );

void           PB_Cainfog2l    ( Int,       Int,       Int,
                                 Int,       Int *,     Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int * );
void           PB_Cinfog2l     ( Int,       Int,       Int *,
                                 Int,       Int,       Int,
                                 Int,       Int *,     Int *,
                                 Int *,     Int * );
Int            PB_Cg2lrem      ( Int,       Int,       Int,
                                 Int,       Int,       Int );
Int            PB_Cindxg2p     ( Int,       Int,       Int,
                                 Int,       Int,       Int );
Int            PB_Cnumroc      ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_Cnpreroc     ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_Cnnxtroc     ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );

void           PB_Cconjg       ( PBTYP_T *, char *,    char * );


void           PB_Cwarn        ( Int,       Int,       char *,
                                 char *,    ... );
void           PB_Cabort       ( Int,       char *,    Int );
void           PB_Cchkmat      ( Int,       char *,    char *,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int *,     Int,       Int * );
void           PB_Cchkvec      ( Int,       char *,    char *,
                                 Int,       Int ,       Int ,
                                 Int,       Int *,     Int,
                                 Int,       Int * );

char *         PB_Cmalloc      ( Int );
char *         PB_Cgetbuf      ( char *,    Int );

PBTYP_T *      PB_Citypeset    ( void );
PBTYP_T *      PB_Cstypeset    ( void );
PBTYP_T *      PB_Cetypeset    ( void );
PBTYP_T *      PB_Cdtypeset    ( void );
PBTYP_T *      PB_Cetypeset    ( void );
PBTYP_T *      PB_Cctypeset    ( void );
PBTYP_T *      PB_Cytypeset    ( void );
PBTYP_T *      PB_Cztypeset    ( void );
PBTYP_T *      PB_Cytypeset    ( void );

Int            pilaenv_        ( Int *,     F_CHAR_T );
char *         PB_Ctop         ( Int *,     char *,    char *,
                                 char * );

void           PB_CVMinit      ( PB_VM_T *, Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int );
Int            PB_CVMnpq       ( PB_VM_T * );
void           PB_CVMcontig    ( PB_VM_T *, Int *,     Int *,
                                 Int *,     Int * );
Int            PB_CVMloc       ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       char *,
                                 char *,    Int );
Int            PB_CVMswp       ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    Int,
                                 char *,    Int,       char *,
                                 Int );
Int            PB_CVMpack      ( PBTYP_T *, PB_VM_T *, char *,
                                 char *,    char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       char *,
                                 char *,    Int );
void           PB_CVMupdate    ( PB_VM_T *, Int,       Int *,
                                 Int * );

void           PB_Cbinfo       ( Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int,       Int,       Int,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int *,     Int *,
                                 Int *,     Int * );

void           PB_Cplaprnt     ( PBTYP_T *, Int,       Int,
                                 char *,    Int,       Int,
                                 Int *,     Int,       Int,
                                 char * );
void           PB_Cplaprn2     ( PBTYP_T *, Int,       Int,
                                 char *,    Int,       Int,
                                 Int *,     Int,       Int,
                                 char *,    Int,       Int );
void           PB_Cprnt        ( char,      Int,       Int,
                                 Int,       char *,    Int,
                                 Int,       char * );

void           PB_Cplapad      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cplapd2      ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    char *,    Int,
                                 Int,       Int * );
void           PB_Cplascal     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cplasca2     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 char *,    Int,       Int,
                                 Int * );
void           PB_Cplacnjg     ( PBTYP_T *, Int,       Int,
                                 char *,    char *,    Int,
                                 Int,       Int * );

void           PB_CInV         ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char * *,  Int *,     Int * );
void           PB_CInV2        ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    Int,       Int * );
void           PB_CInOutV      ( PBTYP_T *, char *,    Int,
                                 Int,       Int *,     Int,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char *,
                                 char * *,  char * *,  Int *,
                                 Int *,     Int *,     Int * );
void           PB_CInOutV2     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int,
                                 Int *,     Int,       char *,
                                 Int,       Int,       Int *,
                                 char *,    char * *,  Int *,
                                 Int *,     Int *,     Int * );
void           PB_COutV        ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       Int *,
                                 Int,       char * *,  Int *,
                                 Int *,     Int * );
void           PB_CGatherV     ( PBTYP_T *, char *,    char *,
                                 Int,       Int,       char *,
                                 Int,       Int,       Int *,
                                 char *,    char * *,  Int *,
                                 Int * );
void           PB_CScatterV    ( PBTYP_T *, char *,    Int,
                                 Int,       char *,    Int,
                                 Int,       Int *,     char *,
                                 char *,    char *,    Int,
                                 Int,       Int *,     char * );
#else

Int            PB_Cgcd         ();
Int            PB_Clcm         ();

void           PB_Cdescset     ();
void           PB_Cdescribe    ();
void           PB_CargFtoC     ();
Int            PB_Cfirstnb     ();
Int            PB_Clastnb      ();
Int            PB_Cspan        ();

void           PB_Cainfog2l    ();
void           PB_Cinfog2l     ();
Int            PB_Cg2lrem      ();
Int            PB_Cindxg2p     ();
Int            PB_Cnumroc      ();
Int            PB_Cnpreroc     ();
Int            PB_Cnnxtroc     ();

void           PB_Cconjg       ();

void           PB_Cwarn        ();
void           PB_Cabort       ();
void           PB_Cchkmat      ();
void           PB_Cchkvec      ();

char *         PB_Cmalloc      ();
char *         PB_Cgetbuf      ();

PBTYP_T *      PB_Citypeset    ();
PBTYP_T *      PB_Cstypeset    ();
PBTYP_T *      PB_Cetypeset    ();
PBTYP_T *      PB_Cdtypeset    ();
PBTYP_T *      PB_Cetypeset    ();
PBTYP_T *      PB_Cctypeset    ();
PBTYP_T *      PB_Cytypeset    ();
PBTYP_T *      PB_Cztypeset    ();
PBTYP_T *      PB_Cytypeset    ();

Int            pilaenv_        ();
char *         PB_Ctop         ();

void           PB_CVMinit      ();
Int            PB_CVMnpq       ();
void           PB_CVMcontig    ();
Int            PB_CVMloc       ();
Int            PB_CVMswp       ();
Int            PB_CVMpack      ();
void           PB_CVMupdate    ();

void           PB_Cbinfo       ();

void           PB_Cplaprnt     ();
void           PB_Cplaprn2     ();
void           PB_Cprnt        ();

void           PB_Cplapad      ();
void           PB_Cplapd2      ();
void           PB_Cplascal     ();
void           PB_Cplasca2     ();
void           PB_Cplacnjg     ();

void           PB_CInV         ();
void           PB_CInV2        ();
void           PB_CInOutV      ();
void           PB_CInOutV2     ();
void           PB_COutV        ();
void           PB_CGatherV     ();
void           PB_CScatterV    ();

#endif
