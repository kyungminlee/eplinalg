      SUBROUTINE YYDOTC( N, DOTC, X, INCX, Y, INCY )
      IMPLICIT NONE
*
*  -- ScaLAPACK tools routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX(KIND=10)            DOTC
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)            X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  YYDOTC is a simple FORTRAN wrapper around the BLAS function
*  YDOTC returning the result in the parameter list instead.
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX(KIND=10)            CTEMP
      INTEGER            I,IX,IY
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
      CTEMP = (0.0d0,0.0d0)
      DOTC = (0.0d0,0.0d0)
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
         DO i = 1,N
            CTEMP = CTEMP + CONJG(X(I))*Y(I)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            CTEMP = CTEMP + CONJG(X(IX))*Y(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DOTC = CTEMP
*
      RETURN
*
*     End of YYDOTC
*
      END
