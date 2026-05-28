*> \brief \b ECABS1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL(KIND=10) FUNCTION ECABS1(Z)
*
*       .. Scalar Arguments ..
*       COMPLEX(KIND=10) Z
*       ..
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ECABS1 computes |Re(.)| + |Im(.)| of a COMPLEX(KIND=10) number
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX(KIND=10)
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
*> \ingroup abs1
*
*  =====================================================================
      REAL(KIND=10) FUNCTION ECABS1(Z)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX(KIND=10) Z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*
      ECABS1 = ABS(REAL(Z, KIND=10)) + ABS(AIMAG(Z))
      RETURN
*
*     End of ECABS1
*
      END
