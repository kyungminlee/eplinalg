*> \brief \b ELAPY3 returns sqrt(x2+y2+z2).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELAPY3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elapy3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elapy3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elapy3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL(KIND=10) FUNCTION ELAPY3( X, Y, Z )
*
*       .. Scalar Arguments ..
*       REAL(KIND=10)   X, Y, Z
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*> unnecessary overflow and unnecessary underflow.
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
*> \endverbatim
*>
*> \param[in] Z
*> \verbatim
*>          Z is REAL(KIND=10)
*>          X, Y and Z specify the values x, y and z.
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
*> \ingroup lapy3
*
*  =====================================================================
      REAL(KIND=10) FUNCTION ELAPY3( X, Y, Z )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   X, Y, Z
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ZERO
      PARAMETER          ( ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      REAL(KIND=10)   W, XABS, YABS, ZABS, HUGEVAL
*     ..
*     .. External Subroutines ..
      REAL(KIND=10)   ELAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      HUGEVAL = ELAMCH( 'Overflow' )
      XABS = ABS( X )
      YABS = ABS( Y )
      ZABS = ABS( Z )
      W = MAX( XABS, YABS, ZABS )
      IF( W.EQ.ZERO .OR. W.GT.HUGEVAL ) THEN
*     W can be zero for max(0,nan,0)
*     adding all three entries together will make sure
*     NaN will not disappear.
         ELAPY3 =  XABS + YABS + ZABS
      ELSE
         ELAPY3 = W*SQRT( ( XABS / W )**2+( YABS / W )**2+             
     +( ZABS / W )**2 )
      END IF
      RETURN
*
*     End of ELAPY3
*
      END
