
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_espmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, EREAL alpha, const EREAL *ap,
		  const EREAL *x, int incx, EREAL beta,
		  EREAL *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_espmv_x, BLAS_ESPMV_X)
 
  (int *uplo, int *n, EREAL *alpha, const EREAL *ap, const EREAL *x,
   int *incx, EREAL *beta, EREAL *y, int *incy, int *prec) {
  BLAS_espmv_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, ap, x,
	       *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
