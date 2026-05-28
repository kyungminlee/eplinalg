
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_esymm_x(enum blas_order_type order, enum blas_side_type side,
		  enum blas_uplo_type uplo, int m, int n,
		  EREAL alpha, const EREAL *a, int lda,
		  const EREAL *b, int ldb, EREAL beta,
		  EREAL *c, int ldc, enum blas_prec_type prec);


extern void FC_FUNC_(blas_esymm_x, BLAS_ESYMM_X)
 
  (int *side, int *uplo, int *m, int *n, EREAL *alpha, const EREAL *a,
   int *lda, const EREAL *b, int *ldb, EREAL *beta, EREAL *c, int *ldc,
   int *prec) {
  BLAS_esymm_x(blas_colmajor, (enum blas_side_type) *side,
	       (enum blas_uplo_type) *uplo, *m, *n, *alpha, a, *lda, b, *ldb,
	       *beta, c, *ldc, (enum blas_prec_type) *prec);
}
