*> \brief \b ELARMM
*
*  Definition:
*  ===========
*
*      REAL(KIND=10) FUNCTION ELARMM( ANORM, BNORM, CNORM )
*
*     .. Scalar Arguments ..
*      REAL(KIND=10)   ANORM, BNORM, CNORM
*     ..
*
*>  \par Purpose:
*  =======
*>
*> \verbatim
*>
*> ELARMM returns a factor s in (0, 1] such that the linear updates
*>
*>    (s * C) - A * (s * B)  and  (s * C) - (s * A) * B
*>
*> cannot overflow, where A, B, and C are matrices of conforming
*> dimensions.
*>
*> This is an auxiliary routine so there is no argument checking.
*> \endverbatim
*
*  Arguments:
*  =========
*
*> \param[in] ANORM
*> \verbatim
*>          ANORM is REAL(KIND=10)
*>          The infinity norm of A. ANORM >= 0.
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] BNORM
*> \verbatim
*>          BNORM is REAL(KIND=10)
*>          The infinity norm of B. BNORM >= 0.
*> \endverbatim
*>
*> \param[in] CNORM
*> \verbatim
*>          CNORM is REAL(KIND=10)
*>          The infinity norm of C. CNORM >= 0.
*> \endverbatim
*>
*>
*  =====================================================================
*>  References:
*>    C. C. Kjelgaard Mikkelsen and L. Karlsson, Blocked Algorithms for
*>    Robust Solution of Triangular Linear Systems. In: International
*>    Conference on Parallel Processing and Applied Mathematics, pages
*>    68--78. Springer, 2017.
*>
*> \ingroup larmm
*  =====================================================================

      REAL(KIND=10) FUNCTION ELARMM( ANORM, BNORM, CNORM )
      IMPLICIT NONE
*     .. Scalar Arguments ..
      REAL(KIND=10)   ANORM, BNORM, CNORM
*     .. Parameters ..
      REAL(KIND=10)   ONE, HALF, FOUR
      PARAMETER          ( ONE = 1.0E0_10, HALF = 0.5E0_10, FOUR = 
     +4.0E0_10 )
*     ..
*     .. Local Scalars ..
       REAL(KIND=10)   BIGNUM, SMLNUM
*     ..
*     .. External Functions ..
      REAL(KIND=10)   ELAMCH
      EXTERNAL           ELAMCH
*     ..
*     .. Executable Statements ..
*
*
*     Determine machine dependent parameters to control overflow.
*
      SMLNUM = ELAMCH( 'Safe minimum' ) / ELAMCH( 'Precision' )
      BIGNUM = ( ONE / SMLNUM ) / FOUR
*
*     Compute a scale factor.
*
      ELARMM = ONE
      IF( BNORM .LE. ONE ) THEN
         IF( ANORM * BNORM .GT. BIGNUM - CNORM ) THEN
            ELARMM = HALF
         END IF
      ELSE
         IF( ANORM .GT. (BIGNUM - CNORM) / BNORM ) THEN
            ELARMM = HALF / BNORM
         END IF
      END IF
      RETURN
*
*     ==== End of ELARMM ====
*
      END
