*> \brief \b ELAPY2 returns sqrt(x2+y2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELAPY2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elapy2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elapy2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elapy2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL(KIND=10) FUNCTION ELAPY2( X, Y )
*
*       .. Scalar Arguments ..
*       REAL(KIND=10)   X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*> overflow and unnecessary underflow.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is REAL(KIND=10)
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is REAL(KIND=10)
*>          X and Y specify the values x and y.
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
*> \ingroup lapy2
*
*  =====================================================================
      REAL(KIND=10) FUNCTION ELAPY2( X, Y )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   X, Y
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ZERO
      PARAMETER          ( ZERO = 0.0E0_10 )
      REAL(KIND=10)   ONE
      PARAMETER          ( ONE = 1.0E0_10 )
*     ..
*     .. Local Scalars ..
      REAL(KIND=10)   W, XABS, YABS, Z, HUGEVAL
      LOGICAL            X_IS_NAN, Y_IS_NAN
*     ..
*     .. External Functions ..
      LOGICAL            EISNAN
      EXTERNAL           EISNAN
*     ..
*     .. External Subroutines ..
      REAL(KIND=10)   ELAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      X_IS_NAN = EISNAN( X )
      Y_IS_NAN = EISNAN( Y )
      IF ( X_IS_NAN ) ELAPY2 = X
      IF ( Y_IS_NAN ) ELAPY2 = Y
      HUGEVAL = ELAMCH( 'Overflow' )
*
      IF ( .NOT.( X_IS_NAN.OR.Y_IS_NAN ) ) THEN
         XABS = ABS( X )
         YABS = ABS( Y )
         W = MAX( XABS, YABS )
         Z = MIN( XABS, YABS )
         IF( Z.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
            ELAPY2 = W
         ELSE
            ELAPY2 = W*SQRT( ONE+( Z / W )**2 )
         END IF
      END IF
      RETURN
*
*     End of ELAPY2
*
      END
