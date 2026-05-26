
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_yhpmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, const void *alpha, const void *ap,
		  const void *x, int incx, const void *beta, void *y,
		  int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_yhpmv_x, BLAS_YHPMV_X)
 
  (int *uplo, int *n, const void *alpha, const void *ap, const void *x,
   int *incx, const void *beta, void *y, int *incy, int *prec) {
  BLAS_yhpmv_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, alpha, ap, x,
	       *incx, beta, y, *incy, (enum blas_prec_type) *prec);
}
