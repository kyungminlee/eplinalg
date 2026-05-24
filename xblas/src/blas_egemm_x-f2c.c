
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_egemm_x(enum blas_order_type order, enum blas_trans_type transa,
		  enum blas_trans_type transb, int m, int n, int k,
		  EREAL alpha, const EREAL *a, int lda, const EREAL *b,
		  int ldb, EREAL beta, EREAL *c, int ldc,
		  enum blas_prec_type prec);


extern void FC_FUNC_(blas_egemm_x, BLAS_EGEMM_X)
 
  (int *transa, int *transb, int *m, int *n, int *k, EREAL *alpha,
   const EREAL *a, int *lda, const EREAL *b, int *ldb, EREAL *beta,
   EREAL *c, int *ldc, int *prec) {
  BLAS_egemm_x(blas_colmajor, (enum blas_trans_type) *transa,
	       (enum blas_trans_type) *transb, *m, *n, *k, *alpha, a, *lda, b,
	       *ldb, *beta, c, *ldc, (enum blas_prec_type) *prec);
}
