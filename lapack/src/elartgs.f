*> \brief \b ELARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELARTGS + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elartgs.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elartgs.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elartgs.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELARTGS( X, Y, SIGMA, CS, SN )
*
*       .. Scalar Arguments ..
*       REAL(KIND=10)        CS, SIGMA, SN, X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELARTGS generates a plane rotation designed to introduce a bulge in
*> Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD
*> problem. X and Y are the top-row entries, and SIGMA is the shift.
*> The computed CS and SN define a plane rotation satisfying
*>
*>    [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],
*>    [ -SN  CS  ]     [    X * Y    ]     [ 0 ]
*>
*> with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the
*> rotation is by PI/2.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is REAL(KIND=10)
*>          The (1,1) entry of an upper bidiagonal matrix.
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is REAL(KIND=10)
*>          The (1,2) entry of an upper bidiagonal matrix.
*> \endverbatim
*>
*> \param[in] SIGMA
*> \verbatim
*>          SIGMA is REAL(KIND=10)
*>          The shift.
*> \endverbatim
*>
*> \param[out] CS
*> \verbatim
*>          CS is REAL(KIND=10)
*>          The cosine of the rotation.
*> \endverbatim
*>
*> \param[out] SN
*> \verbatim
*>          SN is REAL(KIND=10)
*>          The sine of the rotation.
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
*> \ingroup lartgs
*
*  =====================================================================
      SUBROUTINE ELARTGS( X, Y, SIGMA, CS, SN )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)        CS, SIGMA, SN, X, Y
*     ..
*
*  ===================================================================
*
*     .. Parameters ..
      REAL(KIND=10)        NEGONE, ONE, ZERO
      PARAMETER          ( NEGONE = -1.0E0_10, ONE = 1.0E0_10, ZERO = 
     +0.0E0_10 )
*     ..
*     .. Local Scalars ..
      REAL(KIND=10)        R, S, THRESH, W, Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELARTGP
*     ..
*     .. External Functions ..
      REAL(KIND=10)        ELAMCH
      EXTERNAL           ELAMCH
*     .. Executable Statements ..
*
      THRESH = ELAMCH('E')
*
*     Compute the first column of B**T*B - SIGMA^2*I, up to a scale
*     factor.
*
      IF( (SIGMA .EQ. ZERO .AND. ABS(X) .LT. THRESH) .OR.
     $          (ABS(X) .EQ. SIGMA .AND. Y .EQ. ZERO) ) THEN
         Z = ZERO
         W = ZERO
      ELSE IF( SIGMA .EQ. ZERO ) THEN
         IF( X .GE. ZERO ) THEN
            Z = X
            W = Y
         ELSE
            Z = -X
            W = -Y
         END IF
      ELSE IF( ABS(X) .LT. THRESH ) THEN
         Z = -SIGMA*SIGMA
         W = ZERO
      ELSE
         IF( X .GE. ZERO ) THEN
            S = ONE
         ELSE
            S = NEGONE
         END IF
         Z = S * (ABS(X)-SIGMA) * (S+SIGMA/X)
         W = S * Y
      END IF
*
*     Generate the rotation.
*     CALL ELARTGP( Z, W, CS, SN, R ) might seem more natural;
*     reordering the arguments ensures that if Z = 0 then the rotation
*     is by PI/2.
*
      CALL ELARTGP( W, Z, SN, CS, R )
*
      RETURN
*
*     End ELARTGS
*
      END

