*> \brief \b EYSUM1 forms the 1-norm of the complex vector using the true absolute value.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EYSUM1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/eysum1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/eysum1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/eysum1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL(KIND=10) FUNCTION EYSUM1( N, CX, INCX )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       ..
*       .. Array Arguments ..
*       COMPLEX(KIND=10)         CX( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EYSUM1 takes the sum of the absolute values of a complex
*> vector and returns a REAL(KIND=10) result.
*>
*> Based on EYASUM from the Level 1 BLAS.
*> The change is to use the 'genuine' absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements in the vector CX.
*> \endverbatim
*>
*> \param[in] CX
*> \verbatim
*>          CX is COMPLEX(KIND=10) array, dimension (N)
*>          The vector whose elements will be summed.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The spacing between successive values of CX.  INCX > 0.
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
*> \ingroup sum1
*
*> \par Contributors:
*  ==================
*>
*> Nick Higham for use with YLACON.
*
*  =====================================================================
      REAL(KIND=10) FUNCTION EYSUM1( N, CX, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)         CX( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL(KIND=10)   STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      EYSUM1 = 0.0E0_10
      STEMP = 0.0E0_10
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 )
     $   GO TO 20
*
*     CODE FOR INCREMENT NOT EQUAL TO 1
*
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      EYSUM1 = STEMP
      RETURN
*
*     CODE FOR INCREMENT EQUAL TO 1
*
   20 CONTINUE
      DO 30 I = 1, N
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      EYSUM1 = STEMP
      RETURN
*
*     End of EYSUM1
*
      END
