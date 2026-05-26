*> \brief \b EGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by egetc2.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EGESC2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/egesc2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/egesc2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/egesc2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE EGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, N
*       REAL(KIND=10)   SCALE
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), JPIV( * )
*       REAL(KIND=10)   A( LDA, * ), RHS( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EGESC2 solves a system of linear equations
*>
*>           A * X = scale* RHS
*>
*> with a general N-by-N matrix A using the LU factorization with
*> complete pivoting computed by EGETC2.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL(KIND=10) array, dimension (LDA,N)
*>          On entry, the  LU part of the factorization of the n-by-n
*>          matrix A computed by EGETC2:  A = P * L * U * Q
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1, N).
*> \endverbatim
*>
*> \param[in,out] RHS
*> \verbatim
*>          RHS is REAL(KIND=10) array, dimension (N).
*>          On entry, the right hand side vector b.
*>          On exit, the solution vector X.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= i <= N, row i of the
*>          matrix has been interchanged with row IPIV(i).
*> \endverbatim
*>
*> \param[in] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N).
*>          The pivot indices; for 1 <= j <= N, column j of the
*>          matrix has been interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] SCALE
*> \verbatim
*>          SCALE is REAL(KIND=10)
*>          On exit, SCALE contains the scale factor. SCALE is chosen
*>          0 <= SCALE <= 1 to prevent overflow in the solution.
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
*> \ingroup gesc2
*
*> \par Contributors:
*  ==================
*>
*>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*>     Umea University, S-901 87 Umea, Sweden.
*
*  =====================================================================
      SUBROUTINE EGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, N
      REAL(KIND=10)   SCALE
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), JPIV( * )
      REAL(KIND=10)   A( LDA, * ), RHS( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, TWO
      PARAMETER          ( ONE = 1.0E0_10, TWO = 2.0E0_10 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      REAL(KIND=10)   BIGNUM, EPS, SMLNUM, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELASWP, ESCAL
*     ..
*     .. External Functions ..
      INTEGER            IEAMAX
      REAL(KIND=10)   ELAMCH
      EXTERNAL           IEAMAX, ELAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*      Set constant to control overflow
*
      EPS = ELAMCH( 'P' )
      SMLNUM = ELAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Apply permutations IPIV to RHS
*
      CALL ELASWP( 1, RHS, LDA, 1, N-1, IPIV, 1 )
*
*     Solve for L part
*
      DO 20 I = 1, N - 1
         DO 10 J = I + 1, N
            RHS( J ) = RHS( J ) - A( J, I )*RHS( I )
   10    CONTINUE
   20 CONTINUE
*
*     Solve for U part
*
      SCALE = ONE
*
*     Check for scaling
*
      I = IEAMAX( N, RHS, 1 )
      IF( TWO*SMLNUM*ABS( RHS( I ) ).GT.ABS( A( N, N ) ) ) THEN
         TEMP = ( ONE / TWO ) / ABS( RHS( I ) )
         CALL ESCAL( N, TEMP, RHS( 1 ), 1 )
         SCALE = SCALE*TEMP
      END IF
*
      DO 40 I = N, 1, -1
         TEMP = ONE / A( I, I )
         RHS( I ) = RHS( I )*TEMP
         DO 30 J = I + 1, N
            RHS( I ) = RHS( I ) - RHS( J )*( A( I, J )*TEMP )
   30    CONTINUE
   40 CONTINUE
*
*     Apply permutations JPIV to the solution (RHS)
*
      CALL ELASWP( 1, RHS, LDA, 1, N-1, JPIV, -1 )
      RETURN
*
*     End of EGESC2
*
      END
