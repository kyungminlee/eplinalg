      REAL(KIND=10)   FUNCTION PELAMCH( ICTXT, CMACH )
      IMPLICIT NONE
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
      INTEGER            ICTXT
*     ..
*
*  Purpose
*  =======
*
*  PELAMCH determines REAL(KIND=10) machine parameters.
*
*  Arguments
*  =========
*
*  ICTXT   (global input) INTEGER
*          The BLACS context handle in which the computation takes
*          place.
*
*  CMACH   (global input) CHARACTER*1
*          Specifies the value to be returned by PELAMCH:
*          = 'E' or 'e',   PELAMCH := eps
*          = 'S' or 's ,   PELAMCH := sfmin
*          = 'B' or 'b',   PELAMCH := base
*          = 'P' or 'p',   PELAMCH := eps*base
*          = 'N' or 'n',   PELAMCH := t
*          = 'R' or 'r',   PELAMCH := rnd
*          = 'M' or 'm',   PELAMCH := emin
*          = 'U' or 'u',   PELAMCH := rmin
*          = 'L' or 'l',   PELAMCH := emax
*          = 'O' or 'o',   PELAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            IDUMM
      REAL(KIND=10)   TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           EGAMN2D, EGAMX2D
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL(KIND=10)   ELAMCH
      EXTERNAL           ELAMCH, LSAME
*     ..
*     .. Executable Statements ..
*
      TEMP = ELAMCH( CMACH )
      IDUMM = 0
*
      IF( LSAME( CMACH, 'E' ).OR.LSAME( CMACH, 'S' ).OR.
     $    LSAME( CMACH, 'M' ).OR.LSAME( CMACH, 'U' ) ) THEN
         CALL EGAMX2D( ICTXT, 'All', ' ', 1, 1, TEMP, 1, IDUMM,         
     +         IDUMM, -1, -1, IDUMM )
      ELSE IF( LSAME( CMACH, 'L' ).OR.LSAME( CMACH, 'O' ) ) THEN
         CALL EGAMN2D( ICTXT, 'All', ' ', 1, 1, TEMP, 1, IDUMM,         
     +         IDUMM, -1, -1, IDUMM )
      END IF
*
      PELAMCH = TEMP
*
*     End of PELAMCH
*
      END
