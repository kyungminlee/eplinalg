      SUBROUTINE PYLACON( N, V, IV, JV, DESCV, X, IX, JX, DESCX, EST,   
     +                  KASE )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IV, IX, JV, JX, KASE, N
      REAL(KIND=10)               EST
*     ..
*     .. Array Arguments ..
      INTEGER            DESCV( * ), DESCX( * )
      COMPLEX(KIND=10)            V( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PYLACON estimates the 1-norm of a square, complex distributed matrix
*  A. Reverse communication is used for evaluating matrix-vector
*  products. X and V are aligned with the distributed matrix A, this
*  information is implicitly contained within IV, IX, DESCV, and DESCX.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The length of the distributed vectors V and X.  N >= 0.
*
*  V       (local workspace) COMPLEX pointer into the local
*          memory to an array of dimension LOCr(N+MOD(IV-1,MB_V)). On
*          the final return, V = A*W, where EST = norm(V)/norm(W)
*          (W is not returned).
*
*  IV      (global input) INTEGER
*          The row index in the global array V indicating the first
*          row of sub( V ).
*
*  JV      (global input) INTEGER
*          The column index in the global array V indicating the
*          first column of sub( V ).
*
*  DESCV   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix V.
*
*  X       (local input/local output) COMPLEX pointer into the
*          local memory to an array of dimension
*          LOCr(N+MOD(IX-1,MB_X)). On an intermediate return, X
*          should be overwritten by
*                A * X,   if KASE=1,
*                A' * X,  if KASE=2,
*          where A' is the conjugate transpose of A, and PYLACON must
*          be re-called with all the other parameters unchanged.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*
*  EST     (global output) REAL
*          An estimate (a lower bound) for norm(A).
*
*  KASE    (local input/local output) INTEGER
*          On the initial call to PYLACON, KASE should be 0.
*          On an intermediate return, KASE will be 1 or 2, indicating
*          whether X should be overwritten by A * X  or A' * X.
*          On the final return from PYLACON, KASE will again be 0.
*
*  Further Details
*  ===============
*
*  The serial version YLACON has been contributed by Nick Higham,
*  University of Manchester. It was originally named SONEST, dated
*  March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      INTEGER            ITMAX
      PARAMETER          ( ITMAX = 5 )
      REAL(KIND=10)               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0_10, ONE = 1.0E0_10, TWO = 
     +2.0E0_10 )
      COMPLEX(KIND=10)            CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0E0_10, 0.0E0_10 ),              
     +      CONE = ( 1.0E0_10, 0.0E0_10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICTXT, IIVX, IMAXROW, IOFFVX, IROFF, ITER,
     $                   IVXCOL, IVXROW, J, JLAST, JJVX, JUMP, K,
     $                   MYCOL, MYROW, NP, NPCOL, NPROW
      REAL(KIND=10)               ALTSGN, ESTOLD, SAFMIN, TEMP
      COMPLEX(KIND=10)            JLMAX, XMAX
*     ..
*     .. Local Arrays ..
      COMPLEX(KIND=10)            WORK( 2 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, YCOPY, YGEBR2D, YGEBS2D,       
     +             INFOG2L, PYELGET, PYMAX1,                    PEYSUM1,
     + EGEBR2D, EGEBS2D
*     ..
*     .. External Functions ..
      INTEGER            INDXG2L, INDXG2P, INDXL2G, NUMROC
      REAL(KIND=10)               PELAMCH
      EXTERNAL           INDXG2L, INDXG2P, INDXL2G, NUMROC, PELAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
*     Get grid parameters.
*
      EST = ZERO
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      CALL INFOG2L( IX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $              IIVX, JJVX, IVXROW, IVXCOL )
      IF( MYCOL.NE.IVXCOL )
     $   RETURN
      IROFF = MOD( IX-1, DESCX( MB_ ) )
      NP = NUMROC( N+IROFF, DESCX( MB_ ), MYROW, IVXROW, NPROW )
      IF( MYROW.EQ.IVXROW )
     $   NP = NP - IROFF
      IOFFVX = IIVX + (JJVX-1)*DESCX( LLD_ )
*
      SAFMIN = PELAMCH( ICTXT, 'Safe minimum' )
      IF( KASE.EQ.0 ) THEN
         DO 10 I = IOFFVX, IOFFVX+NP-1
            X( I ) = CMPLX( ONE / REAL( N , KIND=10) , KIND=10)
   10    CONTINUE
         KASE = 1
         JUMP = 1
         RETURN
      END IF
*
      GO TO ( 20, 40, 70, 90, 120 )JUMP
*
*     ................ ENTRY   (JUMP = 1)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X
*
   20 CONTINUE
      IF( N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            V( IOFFVX ) = X( IOFFVX )
            EST = ABS( V( IOFFVX ) )
            CALL EGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1 )
         ELSE
            CALL EGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1,       
     +              IVXROW, MYCOL )
         END IF
*        ... QUIT
         GO TO 130
      END IF
      CALL PEYSUM1( N, EST, X, IX, JX, DESCX, 1 )
      IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            CALL EGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1 )
         ELSE
            CALL EGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1,       
     +              IVXROW, MYCOL )
         END IF
      END IF
*
      DO 30 I = IOFFVX, IOFFVX+NP-1
         IF( ABS( X( I ) ).GT.SAFMIN ) THEN
            X( I ) = X( I ) / CMPLX( ABS( X( I ) ) , KIND=10)
         ELSE
            X( I ) = CONE
         END IF
   30 CONTINUE
      KASE = 2
      JUMP = 2
      RETURN
*
*     ................ ENTRY   (JUMP = 2)
*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X
*
   40 CONTINUE
      CALL PYMAX1( N, XMAX, J, X, IX, JX, DESCX, 1 )
      IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            WORK( 1 ) = XMAX
            WORK( 2 ) = CMPLX( REAL( J , KIND=10) , KIND=10)
            CALL YGEBS2D( ICTXT, 'Columnwise', ' ', 2, 1, WORK, 2 )
         ELSE
            CALL YGEBR2D( ICTXT, 'Columnwise', ' ', 2, 1, WORK, 2,      
     +               IVXROW, MYCOL )
            XMAX = WORK( 1 )
            J = NINT( REAL( WORK( 2 ) , KIND=10) )
         END IF
      END IF
      ITER = 2
*
*     MAIN LOOP - ITERATIONS 2, 3,...,ITMAX
*
   50 CONTINUE
      DO 60 I = IOFFVX, IOFFVX+NP-1
         X( I ) = CZERO
   60 CONTINUE
      IMAXROW = INDXG2P( J, DESCX( MB_ ), MYROW, DESCX( RSRC_ ), NPROW )
      IF( MYROW.EQ.IMAXROW ) THEN
         I = INDXG2L( J, DESCX( MB_ ), MYROW, DESCX( RSRC_ ), NPROW )
         X( I ) = CONE
      END IF
      KASE = 1
      JUMP = 3
      RETURN
*
*     ................ ENTRY   (JUMP = 3)
*     X HAS BEEN OVERWRITTEN BY A*X
*
   70 CONTINUE
      CALL YCOPY( NP, X( IOFFVX ), 1, V( IOFFVX ), 1 )
      ESTOLD = EST
      CALL PEYSUM1( N, EST, V, IV, JV, DESCV, 1 )
      IF( DESCV( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            CALL EGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1 )
         ELSE
            CALL EGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, EST, 1,       
     +              IVXROW, MYCOL )
         END IF
      END IF
*
*     TEST FOR CYCLING
      IF( EST.LE.ESTOLD )
     $   GO TO 100
*
      DO 80 I = IOFFVX, IOFFVX+NP-1
         IF( ABS( X( I ) ).GT.SAFMIN ) THEN
            X( I ) = X( I ) / CMPLX( ABS( X( I ) ) , KIND=10)
         ELSE
            X( I ) = CONE
         END IF
   80 CONTINUE
      KASE = 2
      JUMP = 4
      RETURN
*
*     ................ ENTRY   (JUMP = 4)
*     X HAS BEEN OVERWRITTEN BY CTRANS(A)*X
*
   90 CONTINUE
      JLAST = J
      CALL PYMAX1( N, XMAX, J, X, IX, JX, DESCX, 1 )
      IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            WORK( 1 ) = XMAX
            WORK( 2 ) = CMPLX( REAL( J , KIND=10) , KIND=10)
            CALL YGEBS2D( ICTXT, 'Columnwise', ' ', 2, 1, WORK, 2 )
         ELSE
            CALL YGEBR2D( ICTXT, 'Columnwise', ' ', 2, 1, WORK, 2,      
     +               IVXROW, MYCOL )
            XMAX = WORK( 1 )
            J = NINT( REAL( WORK( 2 ) , KIND=10) )
         END IF
      END IF
      CALL PYELGET( 'Columnwise', ' ', JLMAX, X, JLAST, JX, DESCX )
      IF( ( REAL( JLMAX , KIND=10).NE.ABS( REAL( XMAX , KIND=10) ) 
     +).AND.     ( ITER.LT.ITMAX                        ) ) THEN
         ITER = ITER + 1
         GO TO 50
      END IF
*
*     ITERATION COMPLETE.  FINAL STAGE.
*
  100 CONTINUE
      DO 110 I = IOFFVX, IOFFVX+NP-1
         K = INDXL2G( I-IOFFVX+IIVX, DESCX( MB_ ), MYROW,
     $                DESCX( RSRC_ ), NPROW )-IX+1
         IF( MOD( K, 2 ).EQ.0 ) THEN
            ALTSGN = -ONE
         ELSE
            ALTSGN = ONE
         END IF
         X( I ) = CMPLX( ALTSGN*( ONE+REAL( K-1 , KIND=10) / REAL( N-1 ,
     + KIND=10) ) , KIND=10)
  110 CONTINUE
      KASE = 1
      JUMP = 5
      RETURN
*
*     ................ ENTRY   (JUMP = 5)
*     X HAS BEEN OVERWRITTEN BY A*X
*
  120 CONTINUE
      CALL PEYSUM1( N, TEMP, X, IX, JX, DESCX, 1 )
      IF( DESCX( M_ ).EQ.1 .AND. N.EQ.1 ) THEN
         IF( MYROW.EQ.IVXROW ) THEN
            CALL EGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, TEMP, 1 )
         ELSE
            CALL EGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, TEMP, 1,      
     +               IVXROW, MYCOL )
         END IF
      END IF
      TEMP = TWO*( TEMP / REAL( 3*N , KIND=10) )
      IF( TEMP.GT.EST ) THEN
         CALL YCOPY( NP, X( IOFFVX ), 1, V( IOFFVX ), 1 )
         EST = TEMP
      END IF
*
  130 CONTINUE
      KASE = 0
*
      RETURN
*
*     End of PYLACON
*
      END
