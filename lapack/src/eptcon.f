*> \brief \b EPTCON
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EPTCON + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/eptcon.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/eptcon.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/eptcon.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE EPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, N
*       REAL(KIND=10)   ANORM, RCOND
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   D( * ), E( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EPTCON computes the reciprocal of the condition number (in the
*> 1-norm) of a real symmetric positive definite tridiagonal matrix
*> using the factorization A = L*D*L**T or A = U**T*D*U computed by
*> EPTTRF.
*>
*> Norm(inv(A)) is computed by a direct method, and the reciprocal of
*> the condition number is computed as
*>              RCOND = 1 / (ANORM * norm(inv(A))).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL(KIND=10) array, dimension (N)
*>          The n diagonal elements of the diagonal matrix D from the
*>          factorization of A, as computed by EPTTRF.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL(KIND=10) array, dimension (N-1)
*>          The (n-1) off-diagonal elements of the unit bidiagonal factor
*>          U or L from the factorization of A,  as computed by EPTTRF.
*> \endverbatim
*>
*> \param[in] ANORM
*> \verbatim
*>          ANORM is REAL(KIND=10)
*>          The 1-norm of the original matrix A.
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL(KIND=10)
*>          The reciprocal of the condition number of the matrix A,
*>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
*>          1-norm of inv(A) computed in this routine.
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
*> \ingroup ptcon
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The method used is described in Nicholas J. Higham, "Efficient
*>  Algorithms for Computing the Condition Number of a Tridiagonal
*>  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE EPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, N
      REAL(KIND=10)   ANORM, RCOND
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   D( * ), E( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX
      REAL(KIND=10)   AINVNM
*     ..
*     .. External Functions ..
      INTEGER            IEAMAX
      EXTERNAL           IEAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( ANORM.LT.ZERO ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'EPTCON', -INFO )
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
*     Check that D(1:N) is positive.
*
      DO 10 I = 1, N
         IF( D( I ).LE.ZERO )
     $      RETURN
   10 CONTINUE
*
*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
*
*        m(i,j) =  abs(A(i,j)), i = j,
*        m(i,j) = -abs(A(i,j)), i .ne. j,
*
*     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.
*
*     Solve M(L) * x = e.
*
      WORK( 1 ) = ONE
      DO 20 I = 2, N
         WORK( I ) = ONE + WORK( I-1 )*ABS( E( I-1 ) )
   20 CONTINUE
*
*     Solve D * M(L)**T * x = b.
*
      WORK( N ) = WORK( N ) / D( N )
      DO 30 I = N - 1, 1, -1
         WORK( I ) = WORK( I ) / D( I ) + WORK( I+1 )*ABS( E( I ) )
   30 CONTINUE
*
*     Compute AINVNM = max(x(i)), 1<=i<=n.
*
      IX = IEAMAX( N, WORK, 1 )
      AINVNM = ABS( WORK( IX ) )
*
*     Compute the reciprocal condition number.
*
      IF( AINVNM.NE.ZERO )
     $   RCOND = ( ONE / AINVNM ) / ANORM
*
      RETURN
*
*     End of EPTCON
*
      END
