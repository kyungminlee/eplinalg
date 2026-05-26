*> \brief \b YESCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE YESCAL(N,DA,ZX,INCX)
*
*       .. Scalar Arguments ..
*       REAL(KIND=10) DA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       COMPLEX(KIND=10) ZX(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    YESCAL scales a vector by a constant.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in] DA
*> \verbatim
*>          DA is REAL(KIND=10)
*>           On entry, DA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in,out] ZX
*> \verbatim
*>          ZX is COMPLEX(KIND=10) array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of ZX
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
*> \ingroup scal
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, 3/11/78.
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE YESCAL(N,DA,ZX,INCX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10) DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10) ZX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,NINCX
*     .. Parameters ..
      REAL(KIND=10) ONE
      PARAMETER (ONE=1.0E0_10)
*     ..
*     .. Intrinsic Functions ..
*     ..
      IF (N.LE.0 .OR. INCX.LE.0 .OR. DA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DO I = 1,N
            ZX(I) = CMPLX(DA*REAL(ZX(I), KIND=10),DA*AIMAG(ZX(I)), 
     +KIND=10)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            ZX(I) = CMPLX(DA*REAL(ZX(I), KIND=10),DA*AIMAG(ZX(I)), 
     +KIND=10)
         END DO
      END IF
      RETURN
*
*     End of YESCAL
*
      END
