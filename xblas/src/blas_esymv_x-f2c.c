
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_esymv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, EREAL alpha, const EREAL *a, int lda,
		  const EREAL *x, int incx, EREAL beta,
		  EREAL *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_esymv_x, BLAS_ESYMV_X)
 
  (int *uplo, int *n, EREAL *alpha, const EREAL *a, int *lda,
   const EREAL *x, int *incx, EREAL *beta, EREAL *y, int *incy,
   int *prec) {
  BLAS_esymv_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
	       *lda, x, *incx, *beta, y, *incy, (enum blas_prec_type) *prec);
}
