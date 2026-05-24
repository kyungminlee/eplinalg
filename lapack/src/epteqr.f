*> \brief \b EPTEQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download EPTEQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/epteqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/epteqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/epteqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE EPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          COMPZ
*       INTEGER            INFO, LDZ, N
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> EPTEQR computes all eigenvalues and, optionally, eigenvectors of a
*> symmetric positive definite tridiagonal matrix by first factoring the
*> matrix using EPTTRF, and then calling EBDSQR to compute the singular
*> values of the bidiagonal factor.
*>
*> This routine computes the eigenvalues of the positive definite
*> tridiagonal matrix to high relative accuracy.  This means that if the
*> eigenvalues range over many orders of magnitude in size, then the
*> small eigenvalues and corresponding eigenvectors will be computed
*> more accurately than, for example, with the standard QR method.
*>
*> The eigenvectors of a full or band symmetric positive definite matrix
*> can also be found if ESYTRD, ESPTRD, or ESBTRD has been used to
*> reduce this matrix to tridiagonal form. (The reduction to tridiagonal
*> form, however, may preclude the possibility of obtaining high
*> relative accuracy in the small eigenvalues of the original matrix, if
*> these eigenvalues range over many orders of magnitude.)
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N':  Compute eigenvalues only.
*>          = 'V':  Compute eigenvectors of original symmetric
*>                  matrix also.  Array Z contains the orthogonal
*>                  matrix used to reduce the original matrix to
*>                  tridiagonal form.
*>          = 'I':  Compute eigenvectors of tridiagonal matrix also.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is REAL(KIND=10) array, dimension (N)
*>          On entry, the n diagonal elements of the tridiagonal
*>          matrix.
*>          On normal exit, D contains the eigenvalues, in descending
*>          order.
*> \endverbatim
*>
*> \param[in,out] E
*> \verbatim
*>          E is REAL(KIND=10) array, dimension (N-1)
*>          On entry, the (n-1) subdiagonal elements of the tridiagonal
*>          matrix.
*>          On exit, E has been destroyed.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is REAL(KIND=10) array, dimension (LDZ, N)
*>          On entry, if COMPZ = 'V', the orthogonal matrix used in the
*>          reduction to tridiagonal form.
*>          On exit, if COMPZ = 'V', the orthonormal eigenvectors of the
*>          original symmetric matrix;
*>          if COMPZ = 'I', the orthonormal eigenvectors of the
*>          tridiagonal matrix.
*>          If INFO > 0 on exit, Z contains the eigenvectors associated
*>          with only the stored eigenvalues.
*>          If  COMPZ = 'N', then Z is not referenced.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= 1, and if
*>          COMPZ = 'V' or 'I', LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL(KIND=10) array, dimension (4*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          > 0:  if INFO = i, and i is:
*>                <= N  the Cholesky factorization of the matrix could
*>                      not be performed because the leading principal
*>                      minor of order i was not positive.
*>                > N   the SVD algorithm failed to converge;
*>                      if INFO = N+i, i off-diagonal elements of the
*>                      bidiagonal factor did not converge to zero.
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
*> \ingroup pteqr
*
*  =====================================================================
      SUBROUTINE EPTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   D( * ), E( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ZERO, ONE
      PARAMETER          ( ZERO = 0.0E0_10, ONE = 1.0E0_10 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           EBDSQR, ELASET, EPTTRF, XERBLA
*     ..
*     .. Local Arrays ..
      REAL(KIND=10)   C( 1, 1 ), VT( 1, 1 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, NRU
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'EPTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.GT.0 )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
      IF( ICOMPZ.EQ.2 )    CALL ELASET( 'Full', N, N, ZERO, ONE, Z, LDZ 
     +)
*
*     Call EPTTRF to factor the matrix.
*
      CALL EPTTRF( N, D, E, INFO )
      IF( INFO.NE.0 )
     $   RETURN
      DO 10 I = 1, N
         D( I ) = SQRT( D( I ) )
   10 CONTINUE
      DO 20 I = 1, N - 1
         E( I ) = E( I )*D( I )
   20 CONTINUE
*
*     Call EBDSQR to compute the singular values/vectors of the
*     bidiagonal factor.
*
      IF( ICOMPZ.GT.0 ) THEN
         NRU = N
      ELSE
         NRU = 0
      END IF
      CALL EBDSQR( 'Lower', N, 0, NRU, 0, D, E, VT, 1, Z, LDZ, C, 1,    
     +          WORK, INFO )
*
*     Square the singular values.
*
      IF( INFO.EQ.0 ) THEN
         DO 30 I = 1, N
            D( I ) = D( I )*D( I )
   30    CONTINUE
      ELSE
         INFO = N + INFO
      END IF
*
      RETURN
*
*     End of EPTEQR
*
      END
