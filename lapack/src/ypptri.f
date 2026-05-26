*> \brief \b YPPTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YPPTRI + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ypptri.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ypptri.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ypptri.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE YPPTRI( UPLO, N, AP, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       COMPLEX(KIND=10)         AP( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YPPTRI computes the inverse of a complex Hermitian positive definite
*> matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
*> computed by YPPTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangular factor is stored in AP;
*>          = 'L':  Lower triangular factor is stored in AP.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] AP
*> \verbatim
*>          AP is COMPLEX(KIND=10) array, dimension (N*(N+1)/2)
*>          On entry, the triangular factor U or L from the Cholesky
*>          factorization A = U**H*U or A = L*L**H, packed columnwise as
*>          a linear array.  The j-th column of U or L is stored in the
*>          array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
*>
*>          On exit, the upper or lower triangle of the (Hermitian)
*>          inverse of A, overwriting the input factor U or L.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, the (i,i) element of the factor U or L is
*>                zero, and the inverse could not be computed.
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
*> \ingroup pptri
*
*  =====================================================================
      SUBROUTINE YPPTRI( UPLO, N, AP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      COMPLEX(KIND=10)         AP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE
      PARAMETER          ( ONE = 1.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ, JJN
      REAL(KIND=10)   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX(KIND=10)         YDOTC
      EXTERNAL           LSAME, YDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, YESCAL, YHPR, YTPMV, YTPTRI
*     ..
*     .. Intrinsic Functions ..
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'YPPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Invert the triangular Cholesky factor U or L.
*
      CALL YTPTRI( UPLO, 'Non-unit', N, AP, INFO )
      IF( INFO.GT.0 )
     $   RETURN
      IF( UPPER ) THEN
*
*        Compute the product inv(U) * inv(U)**H.
*
         JJ = 0
         DO 10 J = 1, N
            JC = JJ + 1
            JJ = JJ + J
            IF( J.GT.1 )          CALL YHPR( 'Upper', J-1, ONE, AP( JC 
     +), 1, AP )
            AJJ = REAL( AP( JJ ) , KIND=10)
            CALL YESCAL( J, AJJ, AP( JC ), 1 )
   10    CONTINUE
*
      ELSE
*
*        Compute the product inv(L)**H * inv(L).
*
         JJ = 1
         DO 20 J = 1, N
            JJN = JJ + N - J + 1
            AP( JJ ) = REAL( YDOTC( N-J+1, AP( JJ ), 1, AP( JJ ),       
     +    1 ) , KIND=10)
            IF( J.LT.N )          CALL YTPMV( 'Lower', 
     +'Conjugate transpose',                      'Non-unit',           
     +           N-J, AP( JJN ), AP( JJ+1 ), 1 )
            JJ = JJN
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of YPPTRI
*
      END
