*> \brief \b EGECON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EGECON + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/egecon.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/egecon.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/egecon.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE EGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            INFO, LDA, N
*       REAL(KIND=10)   ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       INTEGER            IWORK( * )
*       REAL(KIND=10)   A( LDA, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EGECON estimates the reciprocal of the condition number of a general
*> real matrix A, in either the 1-norm or the infinity-norm, using
*> the LU factorization computed by EGETRF.
*>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the
*> condition number is computed as
*>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies whether the 1-norm condition number or the
*>          infinity-norm condition number is required:
*>          = '1' or 'O':  1-norm;
*>          = 'I':         Infinity-norm.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL(KIND=10) array, dimension (LDA,N)
*>          The factors L and U from the factorization A = P*L*U
*>          as computed by EGETRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] ANORM
*> \verbatim
*>          ANORM is REAL(KIND=10)
*>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
*>          If NORM = 'I', the infinity-norm of the original matrix A.
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL(KIND=10)
*>          The reciprocal of the condition number of the matrix A,
*>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL(KIND=10) array, dimension (4*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>                NaNs are illegal values for ANORM, and they propagate to
*>                the output parameter RCOND.
*>                Infinity is illegal for ANORM, and it propagates to the output
*>                parameter RCOND as 0.
*>          = 1:  if RCOND = NaN, or
*>                   RCOND = Inf, or
*>                   the computed norm of the inverse of A is 0.
*>                In the latter, RCOND = 0 is returned.
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
*> \ingroup gecon
*
*  =====================================================================
      SUBROUTINE EGECON( NORM, N, A, LDA, ANORM, RCOND, WORK, IWORK,    
     +                INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, LDA, N
      REAL(KIND=10)   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      REAL(KIND=10)   A( LDA, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ONENRM
      CHARACTER          NORMIN
      INTEGER            IX, KASE, KASE1
      REAL(KIND=10)   AINVNM, SCALE, SL, SMLNUM, SU, HUGEVAL
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, EISNAN
      INTEGER            IEAMAX
      REAL(KIND=10)   ELAMCH
      EXTERNAL           LSAME, IEAMAX, ELAMCH, EISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELACN2, ELATRS, ERSCL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      HUGEVAL = ELAMCH( 'Overflow' )
*
*     Test the input parameters.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'EGECON', -INFO )
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
      ELSE IF( EISNAN( ANORM ) ) THEN
         RCOND = ANORM
         INFO = -5
         RETURN
      ELSE IF( ANORM.GT.HUGEVAL ) THEN
         INFO = -5
         RETURN
      END IF
*
      SMLNUM = ELAMCH( 'Safe minimum' )
*
*     Estimate the norm of inv(A).
*
      AINVNM = ZERO
      NORMIN = 'N'
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   10 CONTINUE
      CALL ELACN2( N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(L).
*
            CALL ELATRS( 'Lower', 'No transpose', 'Unit', NORMIN, N,    
     +                A,                    LDA, WORK, SL, WORK( 2*N+1 
     +), INFO )
*
*           Multiply by inv(U).
*
            CALL ELATRS( 'Upper', 'No transpose', 'Non-unit', NORMIN,   
     +                 N,                    A, LDA, WORK, SU, WORK( 
     +3*N+1 ), INFO )
         ELSE
*
*           Multiply by inv(U**T).
*
            CALL ELATRS( 'Upper', 'Transpose', 'Non-unit', NORMIN, N,   
     +                 A,                    LDA, WORK, SU, WORK( 3*N+1 
     +), INFO )
*
*           Multiply by inv(L**T).
*
            CALL ELATRS( 'Lower', 'Transpose', 'Unit', NORMIN, N, A,    
     +                LDA, WORK, SL, WORK( 2*N+1 ), INFO )
         END IF
*
*        Divide X by 1/(SL*SU) if doing so will not cause overflow.
*
         SCALE = SL*SU
         NORMIN = 'Y'
         IF( SCALE.NE.ONE ) THEN
            IX = IEAMAX( N, WORK, 1 )
            IF( SCALE.LT.ABS( WORK( IX ) )*SMLNUM .OR. SCALE.EQ.ZERO )
     $         GO TO 20
            CALL ERSCL( N, SCALE, WORK, 1 )
         END IF
         GO TO 10
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO ) THEN
         RCOND = ( ONE / AINVNM ) / ANORM
      ELSE
         INFO = 1
         RETURN
      END IF
*
*     Check for NaNs and Infs
*
      IF( EISNAN( RCOND ) .OR. RCOND.GT.HUGEVAL )    INFO = 1
*
   20 CONTINUE
      RETURN
*
*     End of EGECON
*
      END
