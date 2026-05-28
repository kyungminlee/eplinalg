*> \brief \b YLA_GERPVGRW multiplies a square real matrix by a complex matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YLA_GERPVGRW + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/yla_gerpvgrw.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/yla_gerpvgrw.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/yla_gerpvgrw.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL(KIND=10) FUNCTION YLA_GERPVGRW( N, NCOLS, A, LDA, AF,
*                                               LDAF )
*
*       .. Scalar Arguments ..
*       INTEGER            N, NCOLS, LDA, LDAF
*       ..
*       .. Array Arguments ..
*       COMPLEX(KIND=10)         A( LDA, * ), AF( LDAF, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>
*> YLA_GERPVGRW computes the reciprocal pivot growth factor
*> norm(A)/norm(U). The "max absolute element" norm is used. If this is
*> much less than 1, the stability of the LU factorization of the
*> (equilibrated) matrix A could be poor. This also means that the
*> solution X, estimated condition numbers, and error bounds could be
*> unreliable.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>     The number of linear equations, i.e., the order of the
*>     matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NCOLS
*> \verbatim
*>          NCOLS is INTEGER
*>     The number of columns of the matrix A. NCOLS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX(KIND=10) array, dimension (LDA,N)
*>     On entry, the N-by-N matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>     The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] AF
*> \verbatim
*>          AF is COMPLEX(KIND=10) array, dimension (LDAF,N)
*>     The factors L and U from the factorization
*>     A = P*L*U as computed by YGETRF.
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>     The leading dimension of the array AF.  LDAF >= max(1,N).
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
*> \ingroup la_gerpvgrw
*
*  =====================================================================
      REAL(KIND=10) FUNCTION YLA_GERPVGRW( N, NCOLS, A, LDA, AF,        
     +                                 LDAF )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            N, NCOLS, LDA, LDAF
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)         A( LDA, * ), AF( LDAF, * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
      REAL(KIND=10)   AMAX, UMAX, RPVGRW
      COMPLEX(KIND=10)         ZDUM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, ABS
*     ..
*     .. Statement Functions ..
      REAL(KIND=10)   CABS1
*     ..
*     .. Statement Function Definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM , KIND=10) ) + ABS( AIMAG( ZDUM ) 
     +)
*     ..
*     .. Executable Statements ..
*
      RPVGRW = 1.0E0_10

      DO J = 1, NCOLS
         AMAX = 0.0E0_10
         UMAX = 0.0E0_10
         DO I = 1, N
            AMAX = MAX( CABS1( A( I, J ) ), AMAX )
         END DO
         DO I = 1, J
            UMAX = MAX( CABS1( AF( I, J ) ), UMAX )
         END DO
         IF ( UMAX /= 0.0E0_10 ) THEN
            RPVGRW = MIN( AMAX / UMAX, RPVGRW )
         END IF
      END DO
      YLA_GERPVGRW = RPVGRW
*
*     End of YLA_GERPVGRW
*
      END
