*> \brief \b EPTTRS
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EPTTRS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/epttrs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/epttrs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/epttrs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE EPTTRS( N, NRHS, D, E, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   B( LDB, * ), D( * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EPTTRS solves a tridiagonal system of the form
*>    A * X = B
*> using the L*D*L**T factorization of A computed by EPTTRF.  D is a
*> diagonal matrix specified in the vector D, L is a unit bidiagonal
*> matrix whose subdiagonal is specified in the vector E, and X and B
*> are N by NRHS matrices.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the tridiagonal matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL(KIND=10) array, dimension (N)
*>          The n diagonal elements of the diagonal matrix D from the
*>          L*D*L**T factorization of A.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL(KIND=10) array, dimension (N-1)
*>          The (n-1) subdiagonal elements of the unit bidiagonal factor
*>          L from the L*D*L**T factorization of A.  E can also be regarded
*>          as the superdiagonal of the unit bidiagonal factor U from the
*>          factorization A = U**T*D*U.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL(KIND=10) array, dimension (LDB,NRHS)
*>          On entry, the right hand side vectors B for the system of
*>          linear equations.
*>          On exit, the solution vectors, X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
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
*> \ingroup pttrs
*
*  =====================================================================
      SUBROUTINE EPTTRS( N, NRHS, D, E, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   B( LDB, * ), D( * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            J, JB, NB
*     ..
*     .. External Functions ..
      INTEGER            ILAENV_EP
      EXTERNAL           ILAENV_EP
*     ..
*     .. External Subroutines ..
      EXTERNAL           EPTTS2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments.
*
      INFO = 0
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'EPTTRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Determine the number of right-hand sides to solve at a time.
*
      IF( NRHS.EQ.1 ) THEN
         NB = 1
      ELSE
         NB = MAX( 1, ILAENV_EP( 1, 'EPTTRS', ' ', N, NRHS, -1, -1 ) )
      END IF
*
      IF( NB.GE.NRHS ) THEN
         CALL EPTTS2( N, NRHS, D, E, B, LDB )
      ELSE
         DO 10 J = 1, NRHS, NB
            JB = MIN( NRHS-J+1, NB )
            CALL EPTTS2( N, JB, D, E, B( 1, J ), LDB )
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of EPTTRS
*
      END
