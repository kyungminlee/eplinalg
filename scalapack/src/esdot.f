      SUBROUTINE ESDOT( N, DOT, X, INCX, Y, INCY )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      REAL(KIND=10)               DOT
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)               X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ESDOT is a simple FORTRAN wrapper around the BLAS function
*  EDOT returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL(KIND=10)               EDOT
      EXTERNAL           EDOT
*     ..
*     .. Executable Statements ..
*
      DOT = EDOT( N, X, INCX, Y, INCY )
*
      RETURN
*
*     End of ESDOT
*
      END
