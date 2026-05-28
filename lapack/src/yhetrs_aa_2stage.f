*> \brief \b YHETRS_AA_2STAGE
*
* @generated from SRC/esytrs_aa_2stage.f, fortran d -> c, Mon Oct 30 11:59:02 2017
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download YHETRS_AA_2STAGE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/yhetrs_aa_2stage.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/yhetrs_aa_2stage.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/yhetrs_aa_2stage.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*      SUBROUTINE YHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB, IPIV, 
*                                   IPIV2, B, LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N, NRHS, LDA, LTB, LDB, INFO
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * ), IPIV2( * )
*       COMPLEX(KIND=10)         A( LDA, * ), TB( * ), B( LDB, * )
*       ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> YHETRS_AA_2STAGE solves a system of linear equations A*X = B with a 
*> hermitian matrix A using the factorization A = U**H*T*U or
*> A = L*T*L**H computed by YHETRF_AA_2STAGE.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U**H*T*U;
*>          = 'L':  Lower triangular, form is A = L*T*L**H.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX(KIND=10) array, dimension (LDA,N)
*>          Details of factors computed by YHETRF_AA_2STAGE.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] TB
*> \verbatim
*>          TB is COMPLEX(KIND=10) array, dimension (LTB)
*>          Details of factors computed by YHETRF_AA_2STAGE.
*> \endverbatim
*>
*> \param[in] LTB
*> \verbatim
*>          LTB is INTEGER
*>          The size of the array TB. LTB >= 4*N.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges as computed by
*>          YHETRF_AA_2STAGE.
*> \endverbatim
*>
*> \param[in] IPIV2
*> \verbatim
*>          IPIV2 is INTEGER array, dimension (N)
*>          Details of the interchanges as computed by
*>          YHETRF_AA_2STAGE.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX(KIND=10) array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup hetrs_aa_2stage
*
*  =====================================================================
      SUBROUTINE YHETRS_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB,      
     +                        IPIV, IPIV2, B, LDB, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, NRHS, LDA, LTB, LDB, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IPIV2( * )
      COMPLEX(KIND=10)         A( LDA, * ), TB( * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
      COMPLEX(KIND=10)         ONE
      PARAMETER          ( ONE = ( 1.0E0_10, 0.0E0_10 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            LDTB, NB
      LOGICAL            UPPER
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           YGBTRS, YLASWP, YTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LTB.LT.( 4*N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'YHETRS_AA_2STAGE', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Read NB and compute LDTB
*
      NB = INT( TB( 1 ) )
      LDTB = LTB/N
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U**H*T*U.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B -> B
*
            CALL YLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (U**H \ B) -> B    [ (U**H \P**T * B) ]
*
            CALL YTRSM( 'L', 'U', 'C', 'U', N-NB, NRHS, ONE, A(1,       
     +            NB+1),                  LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (U**H \P**T * B) ]
*
         CALL YGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB,    
     +            INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (U \ B) -> B   [ U \ (T \ (U**H \P**T * B) ) ]
*
            CALL YTRSM( 'L', 'U', 'N', 'U', N-NB, NRHS, ONE, A(1,       
     +            NB+1),                   LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B -> B  [ P * (U \ (T \ (U**H \P**T * B) )) ]
*
            CALL YLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
*
      ELSE
*
*        Solve A*X = B, where A = L*T*L**H.
*
         IF( N.GT.NB ) THEN
*
*           Pivot, P**T * B -> B
*
            CALL YLASWP( NRHS, B, LDB, NB+1, N, IPIV, 1 )
*
*           Compute (L \ B) -> B    [ (L \P**T * B) ]
*
            CALL YTRSM( 'L', 'L', 'N', 'U', N-NB, NRHS, ONE, A(NB+1,    
     +               1),                  LDA, B(NB+1, 1), LDB)
*
         END IF
*
*        Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL YGBTRS( 'N', N, NB, NB, NRHS, TB, LDTB, IPIV2, B, LDB,    
     +            INFO)
         IF( N.GT.NB ) THEN
*
*           Compute (L**H \ B) -> B   [ L**H \ (T \ (L \P**T * B) ) ]
*
            CALL YTRSM( 'L', 'L', 'C', 'U', N-NB, NRHS, ONE, A(NB+1,    
     +               1),                   LDA, B(NB+1, 1), LDB)
*
*           Pivot, P * B -> B  [ P * (L**H \ (T \ (L \P**T * B) )) ]
*
            CALL YLASWP( NRHS, B, LDB, NB+1, N, IPIV, -1 )
*
         END IF
      END IF
*
      RETURN
*
*     End of YHETRS_AA_2STAGE
*
      END
