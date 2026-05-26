*> \brief \b ELARNV returns a vector of random numbers from a uniform or normal distribution.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELARNV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elarnv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elarnv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elarnv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELARNV( IDIST, ISEED, N, X )
*
*       .. Scalar Arguments ..
*       INTEGER            IDIST, N
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 )
*       REAL(KIND=10)   X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELARNV returns a vector of n random real numbers from a uniform or
*> normal distribution.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] IDIST
*> \verbatim
*>          IDIST is INTEGER
*>          Specifies the distribution of the random numbers:
*>          = 1:  uniform (0,1)
*>          = 2:  uniform (-1,1)
*>          = 3:  normal (0,1)
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry, the seed of the random number generator; the array
*>          elements must be between 0 and 4095, and ISEED(4) must be
*>          odd.
*>          On exit, the seed is updated.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of random numbers to be generated.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL(KIND=10) array, dimension (N)
*>          The generated random numbers.
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
*> \ingroup larnv
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  This routine calls the auxiliary routine ELARUV to generate random
*>  real numbers from a uniform (0,1) distribution, in batches of up to
*>  128 using vectorisable code. The Box-Muller method is used to
*>  transform numbers from a uniform to a normal distribution.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ELARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      REAL(KIND=10)   X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, TWO
      PARAMETER          ( ONE = 1.0E0_10, TWO = 2.0E0_10 )
      INTEGER            LV
      PARAMETER          ( LV = 128 )
      REAL(KIND=10)   TWOPI
      PARAMETER  ( TWOPI = 6.28318530717958647692528676655900576839E0_1
     +0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IL, IL2, IV
*     ..
*     .. Local Arrays ..
      REAL(KIND=10)   U( LV )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          COS, LOG, MIN, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           ELARUV
*     ..
*     .. Executable Statements ..
*
      DO 40 IV = 1, N, LV / 2
         IL = MIN( LV / 2, N-IV+1 )
         IF( IDIST.EQ.3 ) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
*
*        Call ELARUV to generate IL2 numbers from a uniform (0,1)
*        distribution (IL2 <= LV)
*
         CALL ELARUV( ISEED, IL2, U )
*
         IF( IDIST.EQ.1 ) THEN
*
*           Copy generated numbers
*
            DO 10 I = 1, IL
               X( IV+I-1 ) = U( I )
   10       CONTINUE
         ELSE IF( IDIST.EQ.2 ) THEN
*
*           Convert generated numbers to uniform (-1,1) distribution
*
            DO 20 I = 1, IL
               X( IV+I-1 ) = TWO*U( I ) - ONE
   20       CONTINUE
         ELSE IF( IDIST.EQ.3 ) THEN
*
*           Convert generated numbers to normal (0,1) distribution
*
            DO 30 I = 1, IL
               X( IV+I-1 ) = SQRT( -TWO*LOG( U( 2*I-1 ) ) )*
     $                       COS( TWOPI*U( 2*I ) )
   30       CONTINUE
         END IF
   40 CONTINUE
      RETURN
*
*     End of ELARNV
*
      END
