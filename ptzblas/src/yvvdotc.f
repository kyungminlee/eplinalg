      SUBROUTINE YVVDOTC( N, DOT, X, INCX, Y, INCY )
      IMPLICIT NONE
*
*  -- PBLAS auxiliary routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      COMPLEX(KIND=10)         DOT
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  YVVDOTC computes the following dot product:
*
*     dot = dot + x**H * y,
*
*  where x and y are n vectors.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          On entry, N  specifies the length of the  vectors x and y.  N
*          must be at least zero.
*
*  DOT     (input/output) COMPLEX(KIND=10)
*          On exit, DOT is updated with the dot product of the vectors x
*          and y.
*
*  X       (input) COMPLEX(KIND=10) array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCX ) ). Before entry,  the incremented
*          array X must contain the vector x.
*
*  INCX    (input) INTEGER
*          On entry, INCX specifies the increment for the elements of X.
*          INCX must not be zero.
*
*  Y       (input) COMPLEX(KIND=10) array of dimension at least
*          ( 1 + ( n - 1 )*abs( INCY ) ). Before entry,  the incremented
*          array Y must contain the vector y.
*
*  INCY    (input) INTEGER
*          On entry, INCY specifies the increment for the elements of Y.
*          INCY must not be zero.
*
*  -- Written on April 1, 1998 by
*     Antoine Petitet, University  of  Tennessee, Knoxville 37996, USA.
*
*  =====================================================================
*
*     .. Local Scalars ..
      COMPLEX(KIND=10)         DOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           YYDOTC
*     ..
*     .. Executable Statements ..
*
      CALL YYDOTC( N, DOTC, X, INCX, Y, INCY )
      DOT = DOT + DOTC
*
      RETURN
*
*     End of YVVDOTC
*
      END
