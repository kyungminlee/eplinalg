*> \brief \b YPBCON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YPBCON + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ypbcon.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ypbcon.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ypbcon.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE YPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
*                          RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, KD, LDAB, N
*       REAL(KIND=10)   ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   RWORK( * )
*       COMPLEX(KIND=10)         AB( LDAB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YPBCON estimates the reciprocal of the condition number (in the
*> 1-norm) of a complex Hermitian positive definite band matrix using
*> the Cholesky factorization A = U**H*U or A = L*L**H computed by
*> YPBTRF.
*>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangular factor stored in AB;
*>          = 'L':  Lower triangular factor stored in AB.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KD
*> \verbatim
*>          KD is INTEGER
*>          The number of superdiagonals of the matrix A if UPLO = 'U',
*>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is COMPLEX(KIND=10) array, dimension (LDAB,N)
*>          The triangular factor U or L from the Cholesky factorization
*>          A = U**H*U or A = L*L**H of the band matrix A, stored in the
*>          first KD+1 rows of the array.  The j-th column of U or L is
*>          stored in the j-th column of the array AB as follows:
*>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
*>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KD+1.
*> \endverbatim
*>
*> \param[in] ANORM
*> \verbatim
*>          ANORM is REAL(KIND=10)
*>          The 1-norm (or infinity-norm) of the Hermitian band matrix A.
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL(KIND=10)
*>          The reciprocal of the condition number of the matrix A,
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
*>          estimate of the 1-norm of inv(A) computed in this routine.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX(KIND=10) array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL(KIND=10) array, dimension (N)
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
*> \ingroup pbcon
*
*  =====================================================================
      SUBROUTINE YPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,     
     +               RWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KD, LDAB, N
      REAL(KIND=10)   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   RWORK( * )
      COMPLEX(KIND=10)         AB( LDAB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      CHARACTER          NORMIN
      INTEGER            IX, KASE
      REAL(KIND=10)   AINVNM, SCALE, SCALEL, SCALEU, SMLNUM
      COMPLEX(KIND=10)         ZDUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IYAMAX
      REAL(KIND=10)   ELAMCH
      EXTERNAL           LSAME, IYAMAX, ELAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, YERSCL, YLACN2, YLATBS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Statement Functions ..
      REAL(KIND=10)   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM , KIND=10) ) + ABS( AIMAG( ZDUM ) 
     +)
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'YPBCON', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      RCOND = ZERO
      IF( N.EQ.0 ) THEN
         RCOND = ONE
         RETURN
      ELSE IF( ANORM.EQ.ZERO ) THEN
         RETURN
      END IF
*
      SMLNUM = ELAMCH( 'Safe minimum' )
*
*     Estimate the 1-norm of the inverse.
*
      KASE = 0
      NORMIN = 'N'
   10 CONTINUE
      CALL YLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( UPPER ) THEN
*
*           Multiply by inv(U**H).
*
            CALL YLATBS( 'Upper', 'Conjugate transpose', 'Non-unit',    
     +                NORMIN, N, KD, AB, LDAB, WORK, SCALEL, RWORK,     
     +               INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(U).
*
            CALL YLATBS( 'Upper', 'No transpose', 'Non-unit', NORMIN,   
     +                 N,                    KD, AB, LDAB, WORK, SCALEU,
     + RWORK, INFO )
         ELSE
*
*           Multiply by inv(L).
*
            CALL YLATBS( 'Lower', 'No transpose', 'Non-unit', NORMIN,   
     +                 N,                    KD, AB, LDAB, WORK, SCALEL,
     + RWORK, INFO )
            NORMIN = 'Y'
*
*           Multiply by inv(L**H).
*
            CALL YLATBS( 'Lower', 'Conjugate transpose', 'Non-unit',    
     +                NORMIN, N, KD, AB, LDAB, WORK, SCALEU, RWORK,     
     +               INFO )
         END IF
*
*        Multiply by 1/SCALE if doing so will not cause overflow.
*
         SCALE = SCALEL*SCALEU
         IF( SCALE.NE.ONE ) THEN
            IX = IYAMAX( N, WORK, 1 )
            IF( SCALE.LT.CABS1( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL YERSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
   20 CONTINUE
*
      RETURN
*
*     End of YPBCON
*
      END
