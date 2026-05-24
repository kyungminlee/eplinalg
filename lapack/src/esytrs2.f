*> \brief \b ESYTRS2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ESYTRS2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/esytrs2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/esytrs2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/esytrs2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ESYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
*                           WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL(KIND=10)   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ESYTRS2 solves a system of linear equations A*X = B with a real
*> symmetric matrix A using the factorization A = U*D*U**T or
*> A = L*D*L**T computed by ESYTRF and converted by ESYCONV.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL(KIND=10) array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by ESYTRF.
*>          Note that A is input / output. This might be counter-intuitive,
*>          and one may think that A is input only. A is input / output. This
*>          is because, at the start of the subroutine, we permute A in a
*>          "better" form and then we permute A back to its original form at
*>          the end.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by ESYTRF.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL(KIND=10) array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL(KIND=10) array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup hetrs2
*
*  =====================================================================
      SUBROUTINE ESYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,          
     +           WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL(KIND=10)   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE
      PARAMETER          ( ONE = 1.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, J, K, KP
      REAL(KIND=10)   AK, AKM1, AKM1K, BK, BKM1, DENOM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ESCAL, ESYCONV, ESWAP, ETRSM,                  
     +  XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ESYTRS2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Convert A
*
      CALL ESYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U*D*U**T.
*
*       P**T * B
        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )          CALL ESWAP( NRHS, B( K, 1 ), LDB, B( 
     +KP, 1 ), LDB )
            K=K-1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( KP.EQ.-IPIV( K-1 ) )          CALL ESWAP( NRHS, B( K-1, 
     +1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         END IF
        END DO
*
*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
*
        CALL ETRSM('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (U \P**T * B) ]
*
         I=N
         DO WHILE ( I .GE. 1 )
            IF( IPIV(I) .GT. 0 ) THEN
              CALL ESCAL( NRHS, ONE / A( I, I ), B( I, 1 ), LDB )
            ELSEIF ( I .GT. 1) THEN
               IF ( IPIV(I-1) .EQ. IPIV(I) ) THEN
                  AKM1K = WORK(I)
                  AKM1 = A( I-1, I-1 ) / AKM1K
                  AK = A( I, I ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  DO 15 J = 1, NRHS
                     BKM1 = B( I-1, J ) / AKM1K
                     BK = B( I, J ) / AKM1K
                     B( I-1, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I, J ) = ( AKM1*BK-BKM1 ) / DENOM
 15              CONTINUE
               I = I - 1
               ENDIF
            ENDIF
            I = I - 1
         END DO
*
*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
*
         CALL ETRSM('L','U','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
*
        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )          CALL ESWAP( NRHS, B( K, 1 ), LDB, B( 
     +KP, 1 ), LDB )
            K=K+1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K .LT. N .AND. KP.EQ.-IPIV( K+1 ) )          CALL 
     +ESWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L**T.
*
*       P**T * B
        K=1
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )          CALL ESWAP( NRHS, B( K, 1 ), LDB, B( 
     +KP, 1 ), LDB )
            K=K+1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K and -IPIV(K+1).
            KP = -IPIV( K+1 )
            IF( KP.EQ.-IPIV( K ) )          CALL ESWAP( NRHS, B( K+1, 1 
     +), LDB, B( KP, 1 ), LDB )
            K=K+2
         ENDIF
        END DO
*
*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
*
        CALL ETRSM('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (L \P**T * B) ]
*
         I=1
         DO WHILE ( I .LE. N )
            IF( IPIV(I) .GT. 0 ) THEN
              CALL ESCAL( NRHS, ONE / A( I, I ), B( I, 1 ), LDB )
            ELSE
                  AKM1K = WORK(I)
                  AKM1 = A( I, I ) / AKM1K
                  AK = A( I+1, I+1 ) / AKM1K
                  DENOM = AKM1*AK - ONE
                  DO 25 J = 1, NRHS
                     BKM1 = B( I, J ) / AKM1K
                     BK = B( I+1, J ) / AKM1K
                     B( I, J ) = ( AK*BKM1-BK ) / DENOM
                     B( I+1, J ) = ( AKM1*BK-BKM1 ) / DENOM
 25              CONTINUE
                  I = I + 1
            ENDIF
            I = I + 1
         END DO
*
*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
*
        CALL ETRSM('L','L','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
*
        K=N
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           1 x 1 diagonal block
*           Interchange rows K and IPIV(K).
            KP = IPIV( K )
            IF( KP.NE.K )          CALL ESWAP( NRHS, B( K, 1 ), LDB, B( 
     +KP, 1 ), LDB )
            K=K-1
         ELSE
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K).
            KP = -IPIV( K )
            IF( K.GT.1 .AND. KP.EQ.-IPIV( K-1 ) )          CALL ESWAP( 
     +NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            K=K-2
         ENDIF
        END DO
*
      END IF
*
*     Revert A
*
      CALL ESYCONV( UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO )
*
      RETURN
*
*     End of ESYTRS2
*
      END
