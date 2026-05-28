
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ytbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, const void *alpha, const void *t, int ldt,
		  void *x, int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ytbsv_x, BLAS_YTBSV_X)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, const void *alpha,
   const void *t, int *ldt, void *x, int *incx, int *prec) {
  BLAS_ytbsv_x(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *k, alpha, t, *ldt, x, *incx, (enum blas_prec_type) *prec);
}
