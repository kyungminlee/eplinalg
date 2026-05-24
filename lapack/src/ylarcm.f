*> \brief \b YLARCM copies all or part of a real two-dimensional array to a complex array.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YLARCM + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ylarcm.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ylarcm.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ylarcm.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE YLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, LDB, LDC, M, N
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   A( LDA, * ), RWORK( * )
*       COMPLEX(KIND=10)         B( LDB, * ), C( LDC, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YLARCM performs a very simple matrix-matrix multiplication:
*>          C := A * B,
*> where A is M by M and real; B is M by N and complex;
*> C is M by N and complex.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A and of the matrix C.
*>          M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns and rows of the matrix B and
*>          the number of columns of the matrix C.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is REAL(KIND=10) array, dimension (LDA, M)
*>          On entry, A contains the M by M matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >=max(1,M).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX(KIND=10) array, dimension (LDB, N)
*>          On entry, B contains the M by N matrix B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >=max(1,M).
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is COMPLEX(KIND=10) array, dimension (LDC, N)
*>          On exit, C contains the M by N matrix C.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >=max(1,M).
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL(KIND=10) array, dimension (2*M*N)
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
*> \ingroup larcm
*
*  =====================================================================
      SUBROUTINE YLARCM( M, N, A, LDA, B, LDB, C, LDC, RWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDB, LDC, M, N
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   A( LDA, * ), RWORK( * )
      COMPLEX(KIND=10)         B( LDB, * ), C( LDC, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. External Subroutines ..
      EXTERNAL           EGEMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
*
      DO 20 J = 1, N
         DO 10 I = 1, M
            RWORK( ( J-1 )*M+I ) = REAL( B( I, J ) , KIND=10)
   10    CONTINUE
   20 CONTINUE
*
      L = M*N + 1
      CALL EGEMM( 'N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO,       
     +      RWORK( L ), M )
      DO 40 J = 1, N
         DO 30 I = 1, M
            C( I, J ) = RWORK( L+( J-1 )*M+I-1 )
   30    CONTINUE
   40 CONTINUE
*
      DO 60 J = 1, N
         DO 50 I = 1, M
            RWORK( ( J-1 )*M+I ) = AIMAG( B( I, J ) )
   50    CONTINUE
   60 CONTINUE
      CALL EGEMM( 'N', 'N', M, N, M, ONE, A, LDA, RWORK, M, ZERO,       
     +      RWORK( L ), M )
      DO 80 J = 1, N
         DO 70 I = 1, M
            C( I, J ) = CMPLX( REAL( C( I, J ) , KIND=10),              
     +     RWORK( L+( J-1 )*M+I-1 ) , KIND=10)
   70    CONTINUE
   80 CONTINUE
*
      RETURN
*
*     End of YLARCM
*
      END
