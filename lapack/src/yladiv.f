*> \brief \b YLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YLADIV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/yladiv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/yladiv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/yladiv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       COMPLEX(KIND=10)     FUNCTION YLADIV( X, Y )
*
*       .. Scalar Arguments ..
*       COMPLEX(KIND=10)         X, Y
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*> will not overflow on an intermediary step unless the results
*> overflows.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] X
*> \verbatim
*>          X is COMPLEX(KIND=10)
*> \endverbatim
*>
*> \param[in] Y
*> \verbatim
*>          Y is COMPLEX(KIND=10)
*>          The complex scalars X and Y.
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
      COMPLEX(KIND=10)     FUNCTION YLADIV( X, Y )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX(KIND=10)         X, Y
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL(KIND=10)   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELADIV
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
      CALL ELADIV( REAL( X , KIND=10), AIMAG( X ), REAL( Y , KIND=10), 
     +AIMAG( Y ), ZR,              ZI )
      YLADIV = CMPLX( ZR, ZI , KIND=10)
*
      RETURN
*
*     End of YLADIV
*
      END
