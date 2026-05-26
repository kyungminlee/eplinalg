*> \brief \b ELABAD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELABAD + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elabad.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elabad.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elabad.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELABAD( SMALL, LARGE )
*
*       .. Scalar Arguments ..
*       REAL(KIND=10)   LARGE, SMALL
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELABAD is a no-op and kept for compatibility reasons. It used
*> to correct the overflow/underflow behavior of machines that
*> are not IEEE-754 compliant.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in,out] SMALL
*> \verbatim
*>          SMALL is REAL(KIND=10)
*>          On entry, the underflow threshold as computed by ELAMCH.
*>          On exit, the unchanged value SMALL.
*> \endverbatim
*>
*> \param[in,out] LARGE
*> \verbatim
*>          LARGE is REAL(KIND=10)
*>          On entry, the overflow threshold as computed by ELAMCH.
*>          On exit, the unchanged value LARGE.
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
*> \ingroup labad
*
*  =====================================================================
      SUBROUTINE ELABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   LARGE, SMALL
*     ..
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
*      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
*         SMALL = SQRT( SMALL )
*         LARGE = SQRT( LARGE )
*      END IF
*
      RETURN
*
*     End of ELABAD
*
      END
