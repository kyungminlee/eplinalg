*> \brief \b YGTCON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YGTCON + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ygtcon.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ygtcon.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ygtcon.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE YGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,
*                          WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            INFO, N
*       REAL(KIND=10)   ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX(KIND=10)         D( * ), DL( * ), DU( * ), DU2( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YGTCON estimates the reciprocal of the condition number of a complex
*> tridiagonal matrix A using the LU factorization as computed by
*> YGTTRF.
*>
*> An estimate is obtained for norm(inv(A)), and the reciprocal of the
*> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
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
*> \param[in] DL
*> \verbatim
*>          DL is COMPLEX(KIND=10) array, dimension (N-1)
*>          The (n-1) multipliers that define the matrix L from the
*>          LU factorization of A as computed by YGTTRF.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is COMPLEX(KIND=10) array, dimension (N)
*>          The n diagonal elements of the upper triangular matrix U from
*>          the LU factorization of A.
*> \endverbatim
*>
*> \param[in] DU
*> \verbatim
*>          DU is COMPLEX(KIND=10) array, dimension (N-1)
*>          The (n-1) elements of the first superdiagonal of U.
*> \endverbatim
*>
*> \param[in] DU2
*> \verbatim
*>          DU2 is COMPLEX(KIND=10) array, dimension (N-2)
*>          The (n-2) elements of the second superdiagonal of U.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          The pivot indices; for 1 <= i <= n, row i of the matrix was
*>          interchanged with row IPIV(i).  IPIV(i) will always be either
*>          i or i+1; IPIV(i) = i indicates a row interchange was not
*>          required.
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
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
*>          estimate of the 1-norm of inv(A) computed in this routine.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX(KIND=10) array, dimension (2*N)
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
*> \ingroup gtcon
*
*  =====================================================================
      SUBROUTINE YGTCON( NORM, N, DL, D, DU, DU2, IPIV, ANORM, RCOND,   
     +                 WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            INFO, N
      REAL(KIND=10)   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX(KIND=10)         D( * ), DL( * ), DU( * ), DU2( * ), 
     +WORK( * )
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
      INTEGER            I, KASE, KASE1
      REAL(KIND=10)   AINVNM
*     ..
*     .. Local Arrays ..
      INTEGER            ISAVE( 3 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, YGTTRS, YLACN2
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      ONENRM = NORM.EQ.'1' .OR. LSAME( NORM, 'O' )
      IF( .NOT.ONENRM .AND. .NOT.LSAME( NORM, 'I' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'YGTCON', -INFO )
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
*     Check that D(1:N) is non-zero.
*
      DO 10 I = 1, N
         IF( D( I ).EQ.CMPLX( ZERO , KIND=10) )       RETURN
   10 CONTINUE
*
      AINVNM = ZERO
      IF( ONENRM ) THEN
         KASE1 = 1
      ELSE
         KASE1 = 2
      END IF
      KASE = 0
   20 CONTINUE
      CALL YLACN2( N, WORK( N+1 ), WORK, AINVNM, KASE, ISAVE )
      IF( KASE.NE.0 ) THEN
         IF( KASE.EQ.KASE1 ) THEN
*
*           Multiply by inv(U)*inv(L).
*
            CALL YGTTRS( 'No transpose', N, 1, DL, D, DU, DU2, IPIV,    
     +                WORK, N, INFO )
         ELSE
*
*           Multiply by inv(L**H)*inv(U**H).
*
            CALL YGTTRS( 'Conjugate transpose', N, 1, DL, D, DU, DU2,   
     +                 IPIV, WORK, N, INFO )
         END IF
         GO TO 20
      END IF
*
*     Compute the estimate of the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
*     End of YGTCON
*
      END
