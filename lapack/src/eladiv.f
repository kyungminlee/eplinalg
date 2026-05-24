*> \brief \b ELADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELADIV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/eladiv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/eladiv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/eladiv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELADIV( A, B, C, D, P, Q )
*
*       .. Scalar Arguments ..
*       REAL(KIND=10)   A, B, C, D, P, Q
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELADIV performs complex division in  real arithmetic
*>
*>                       a + i*b
*>            p + i*q = ---------
*>                       c + i*d
*>
*> The algorithm is due to Michael Baudin and Robert L. Smith
*> and can be found in the paper
*> "A Robust Complex Division in Scilab"
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is REAL(KIND=10)
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is REAL(KIND=10)
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL(KIND=10)
*> \endverbatim
*>
*> \param[in] D
*> \verbatim
*>          D is REAL(KIND=10)
*>          The scalars a, b, c, and d in the above expression.
*> \endverbatim
*>
*> \param[out] P
*> \verbatim
*>          P is REAL(KIND=10)
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is REAL(KIND=10)
*>          The scalars p and q in the above expression.
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
*> \ingroup ladiv
*
*  =====================================================================
      SUBROUTINE ELADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   BS
      PARAMETER          ( BS = 2.0E0_10 )
      REAL(KIND=10)   HALF
      PARAMETER          ( HALF = 0.5E0_10 )
      REAL(KIND=10)   TWO
      PARAMETER          ( TWO = 2.0E0_10 )
*
*     .. Local Scalars ..
      REAL(KIND=10)   AA, BB, CC, DD, AB, CD, S, OV, UN, BE, EPS
*     ..
*     .. External Functions ..
      REAL(KIND=10)   ELAMCH
      EXTERNAL           ELAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELADIV1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
      AA = A
      BB = B
      CC = C
      DD = D
      AB = MAX( ABS(A), ABS(B) )
      CD = MAX( ABS(C), ABS(D) )
      S = 1.0E0_10

      OV = ELAMCH( 'Overflow threshold' )
      UN = ELAMCH( 'Safe minimum' )
      EPS = ELAMCH( 'Epsilon' )
      BE = BS / (EPS*EPS)

      IF( AB >= HALF*OV ) THEN
         AA = HALF * AA
         BB = HALF * BB
         S  = TWO * S
      END IF
      IF( CD >= HALF*OV ) THEN
         CC = HALF * CC
         DD = HALF * DD
         S  = HALF * S
      END IF
      IF( AB <= UN*BS/EPS ) THEN
         AA = AA * BE
         BB = BB * BE
         S  = S / BE
      END IF
      IF( CD <= UN*BS/EPS ) THEN
         CC = CC * BE
         DD = DD * BE
         S  = S * BE
      END IF
      IF( ABS( D ).LE.ABS( C ) ) THEN
         CALL ELADIV1(AA, BB, CC, DD, P, Q)
      ELSE
         CALL ELADIV1(BB, AA, DD, CC, P, Q)
         Q = -Q
      END IF
      P = P * S
      Q = Q * S
*
      RETURN
*
*     End of ELADIV
*
      END

*> \ingroup ladiv


      SUBROUTINE ELADIV1( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   A, B, C, D, P, Q
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE
      PARAMETER          ( ONE = 1.0E0_10 )
*
*     .. Local Scalars ..
      REAL(KIND=10)   R, T
*     ..
*     .. External Functions ..
      REAL(KIND=10)   ELADIV2
      EXTERNAL           ELADIV2
*     ..
*     .. Executable Statements ..
*
      R = D / C
      T = ONE / (C + D * R)
      P = ELADIV2(A, B, C, D, R, T)
      A = -A
      Q = ELADIV2(B, A, C, D, R, T)
*
      RETURN
*
*     End of ELADIV1
*
      END

*> \ingroup ladiv

      REAL(KIND=10) FUNCTION ELADIV2( A, B, C, D, R, T )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   A, B, C, D, R, T
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ZERO
      PARAMETER          ( ZERO = 0.0E0_10 )
*
*     .. Local Scalars ..
      REAL(KIND=10)   BR
*     ..
*     .. Executable Statements ..
*
      IF( R.NE.ZERO ) THEN
         BR = B * R
         IF( BR.NE.ZERO ) THEN
            ELADIV2 = (A + BR) * T
         ELSE
            ELADIV2 = A * T + (B * T) * R
         END IF
      ELSE
         ELADIV2 = (A + D * (B / C)) * T
      END IF
*
      RETURN
*
*     End of ELADIV2
*
      END
