*> \brief \b YTRSYL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YTRSYL + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ytrsyl.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ytrsyl.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ytrsyl.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE YTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
*                          LDC, SCALE, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          TRANA, TRANB
*       INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
*       REAL(KIND=10)   SCALE
*       ..
*       .. Array Arguments ..
*       COMPLEX(KIND=10)         A( LDA, * ), B( LDB, * ), C( LDC, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YTRSYL solves the complex Sylvester matrix equation:
*>
*>    op(A)*X + X*op(B) = scale*C or
*>    op(A)*X - X*op(B) = scale*C,
*>
*> where op(A) = A or A**H, and A and B are both upper triangular. A is
*> M-by-M and B is N-by-N; the right hand side C and the solution X are
*> M-by-N; and scale is an output scale factor, set <= 1 to avoid
*> overflow in X.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANA
*> \verbatim
*>          TRANA is CHARACTER*1
*>          Specifies the option op(A):
*>          = 'N': op(A) = A    (No transpose)
*>          = 'C': op(A) = A**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] TRANB
*> \verbatim
*>          TRANB is CHARACTER*1
*>          Specifies the option op(B):
*>          = 'N': op(B) = B    (No transpose)
*>          = 'C': op(B) = B**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] ISGN
*> \verbatim
*>          ISGN is INTEGER
*>          Specifies the sign in the equation:
*>          = +1: solve op(A)*X + X*op(B) = scale*C
*>          = -1: solve op(A)*X - X*op(B) = scale*C
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The order of the matrix A, and the number of rows in the
*>          matrices X and C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix B, and the number of columns in the
*>          matrices X and C. N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX(KIND=10) array, dimension (LDA,M)
*>          The upper triangular matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX(KIND=10) array, dimension (LDB,N)
*>          The upper triangular matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX(KIND=10) array, dimension (LDC,N)
*>          On entry, the M-by-N right hand side matrix C.
*>          On exit, C is overwritten by the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M)
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is REAL(KIND=10)
*>          The scale factor, scale, set <= 1 to avoid overflow in X.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          = 1: A and B have common or very close eigenvalues; perturbed
*>               values were used to solve the equation (but the matrices
*>               A and B are unchanged).
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup trsyl
*
*  =====================================================================
      SUBROUTINE YTRSYL( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,   
     +                 LDC, SCALE, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      REAL(KIND=10)   SCALE
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)         A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE
      PARAMETER          ( ONE = 1.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            J, K, L
      REAL(KIND=10)   BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,         
     +           SMLNUM
      COMPLEX(KIND=10)         A11, SUML, SUMR, VEC, X11
*     ..
*     .. Local Arrays ..
      REAL(KIND=10)   DUM( 1 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL(KIND=10)   ELAMCH, YLANGE
      COMPLEX(KIND=10)         YDOTC, YDOTU, YLADIV
      EXTERNAL           LSAME, ELAMCH, YLANGE, YDOTC, YDOTU,           
     +         YLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, YESCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Decode and Test input parameters
*
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
*
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'YTRSYL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Set constants to control overflow
*
      EPS = ELAMCH( 'P' )
      SMLNUM = ELAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      SMLNUM = SMLNUM*REAL( M*N , KIND=10) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*YLANGE( 'M', M, M, A, LDA, DUM ),        
     +EPS*YLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN
*
      IF( NOTRNA .AND. NOTRNB ) THEN
*
*        Solve    A*X + ISGN*X*B = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        bottom-left corner column by column by
*
*            A(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*
*        Where
*                    M                        L-1
*          R(K,L) = SUM [A(K,I)*X(I,L)] +ISGN*SUM [X(K,J)*B(J,L)].
*                  I=K+1                      J=1
*
         DO 30 L = 1, N
            DO 20 K = M, 1, -1
*
               SUML = YDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,           
     +      C( MIN( K+1, M ), L ), 1 )
               SUMR = YDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
*
               SCALOC = ONE
               A11 = A( K, K ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 , KIND=10) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC , KIND=10) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = YLADIV( VEC*CMPLX( SCALOC , KIND=10), A11 )
*
               IF( SCALOC.NE.ONE ) THEN
                  DO 10 J = 1, N
                     CALL YESCAL( M, SCALOC, C( 1, J ), 1 )
   10             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
*
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
*
*        Solve    A**H *X + ISGN*X*B = scale*C.
*
*        The (K,L)th block of X is determined starting from
*        upper-left corner column by column by
*
*            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B(L,L) = C(K,L) - R(K,L)
*
*        Where
*                   K-1                           L-1
*          R(K,L) = SUM [A**H(I,K)*X(I,L)] + ISGN*SUM [X(K,J)*B(J,L)]
*                   I=1                           J=1
*
         DO 60 L = 1, N
            DO 50 K = 1, M
*
               SUML = YDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = YDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
*
               SCALOC = ONE
               A11 = CONJG( A( K, K ) ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 , KIND=10) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC , KIND=10) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
*
               X11 = YLADIV( VEC*CMPLX( SCALOC , KIND=10), A11 )
*
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 J = 1, N
                     CALL YESCAL( M, SCALOC, C( 1, J ), 1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
*
   50       CONTINUE
   60    CONTINUE
*
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
*
*        Solve    A**H*X + ISGN*X*B**H = C.
*
*        The (K,L)th block of X is determined starting from
*        upper-right corner column by column by
*
*            A**H(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
*
*        Where
*                    K-1
*           R(K,L) = SUM [A**H(I,K)*X(I,L)] +
*                    I=1
*                           N
*                     ISGN*SUM [X(K,J)*B**H(L,J)].
*                          J=L+1
*
         DO 90 L = N, 1, -1
            DO 80 K = 1, M
*
               SUML = YDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = YDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,           
     +      B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
*
               SCALOC = ONE
               A11 = CONJG( A( K, K )+SGN*B( L, L ) )
               DA11 = ABS( REAL( A11 , KIND=10) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC , KIND=10) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
*
               X11 = YLADIV( VEC*CMPLX( SCALOC , KIND=10), A11 )
*
               IF( SCALOC.NE.ONE ) THEN
                  DO 70 J = 1, N
                     CALL YESCAL( M, SCALOC, C( 1, J ), 1 )
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
*
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
*
*        Solve    A*X + ISGN*X*B**H = C.
*
*        The (K,L)th block of X is determined starting from
*        bottom-left corner column by column by
*
*           A(K,K)*X(K,L) + ISGN*X(K,L)*B**H(L,L) = C(K,L) - R(K,L)
*
*        Where
*                    M                          N
*          R(K,L) = SUM [A(K,I)*X(I,L)] + ISGN*SUM [X(K,J)*B**H(L,J)]
*                  I=K+1                      J=L+1
*
         DO 120 L = N, 1, -1
            DO 110 K = M, 1, -1
*
               SUML = YDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,           
     +      C( MIN( K+1, M ), L ), 1 )
               SUMR = YDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,           
     +      B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
*
               SCALOC = ONE
               A11 = A( K, K ) + SGN*CONJG( B( L, L ) )
               DA11 = ABS( REAL( A11 , KIND=10) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC , KIND=10) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
*
               X11 = YLADIV( VEC*CMPLX( SCALOC , KIND=10), A11 )
*
               IF( SCALOC.NE.ONE ) THEN
                  DO 100 J = 1, N
                     CALL YESCAL( M, SCALOC, C( 1, J ), 1 )
  100             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
*
  110       CONTINUE
  120    CONTINUE
*
      END IF
*
      RETURN
*
*     End of YTRSYL
*
      END
