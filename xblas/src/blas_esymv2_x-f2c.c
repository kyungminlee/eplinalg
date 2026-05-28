
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_esymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, EREAL alpha, const EREAL *a, int lda,
		   const EREAL *x_head, const EREAL *x_tail, int incx,
		   EREAL beta, EREAL *y, int incy,
		   enum blas_prec_type prec);


extern void FC_FUNC_(blas_esymv2_x, BLAS_ESYMV2_X)
 
  (int *uplo, int *n, EREAL *alpha, const EREAL *a, int *lda,
   const EREAL *x_head, const EREAL *x_tail, int *incx, EREAL *beta,
   EREAL *y, int *incy, int *prec) {
  BLAS_esymv2_x(blas_colmajor, (enum blas_uplo_type) *uplo, *n, *alpha, a,
		*lda, x_head, x_tail, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
