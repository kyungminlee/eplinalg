*> \brief \b ELARFY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INCV, LDC, N
*       REAL(KIND=10)   TAU
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELARFY applies an elementary reflector, or Householder matrix, H,
*> to an n x n symmetric matrix C, from both the left and the right.
*>
*> H is represented in the form
*>
*>    H = I - tau * v * v'
*>
*> where  tau  is a scalar and  v  is a vector.
*>
*> If  tau  is  zero, then  H  is taken to be the unit matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          symmetric matrix C is stored.
*>          = 'U':  Upper triangle
*>          = 'L':  Lower triangle
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix C.  N >= 0.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is REAL(KIND=10) array, dimension
*>                  (1 + (N-1)*abs(INCV))
*>          The vector v as described above.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between successive elements of v.  INCV must
*>          not be zero.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is REAL(KIND=10)
*>          The value tau as described above.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is REAL(KIND=10) array, dimension (LDC, N)
*>          On entry, the matrix C.
*>          On exit, C is overwritten by H * C * H'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C.  LDC >= max( 1, N ).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL(KIND=10) array, dimension (N)
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
*> \ingroup larfy
*
*  =====================================================================
      SUBROUTINE ELARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INCV, LDC, N
      REAL(KIND=10)   TAU
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10, HALF = 
     +0.5E0_10 )
*     ..
*     .. Local Scalars ..
      REAL(KIND=10)   ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           EAXPY, ESYMV, ESYR2
*     ..
*     .. External Functions ..
      REAL(KIND=10)   EDOT
      EXTERNAL           EDOT
*     ..
*     .. Executable Statements ..
*
      IF( TAU.EQ.ZERO )
     $   RETURN
*
*     Form  w:= C * v
*
      CALL ESYMV( UPLO, N, ONE, C, LDC, V, INCV, ZERO, WORK, 1 )
*
      ALPHA = -HALF*TAU*EDOT( N, WORK, 1, V, INCV )
      CALL EAXPY( N, ALPHA, V, INCV, WORK, 1 )
*
*     C := C - v * w' - w * v'
*
      CALL ESYR2( UPLO, N, -TAU, V, INCV, WORK, 1, C, LDC )
*
      RETURN
*
*     End of ELARFY
*
      END
