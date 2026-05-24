
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_etbsv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  enum blas_trans_type trans, enum blas_diag_type diag,
		  int n, int k, EREAL alpha, const EREAL *t, int ldt,
		  EREAL *x, int incx, enum blas_prec_type prec);


extern void FC_FUNC_(blas_etbsv_x, BLAS_ETBSV_X)
 
  (int *uplo, int *trans, int *diag, int *n, int *k, EREAL *alpha,
   const EREAL *t, int *ldt, EREAL *x, int *incx, int *prec) {
  BLAS_etbsv_x(blas_colmajor, (enum blas_uplo_type) *uplo,
	       (enum blas_trans_type) *trans, (enum blas_diag_type) *diag, *n,
	       *k, *alpha, t, *ldt, x, *incx, (enum blas_prec_type) *prec);
}
