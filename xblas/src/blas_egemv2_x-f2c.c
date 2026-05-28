
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_egemv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, EREAL alpha, const EREAL *a, int lda,
		   const EREAL *head_x, const EREAL *tail_x, int incx,
		   EREAL beta, EREAL *y, int incy,
		   enum blas_prec_type prec);


extern void FC_FUNC_(blas_egemv2_x, BLAS_EGEMV2_X)
 
  (int *trans, int *m, int *n, EREAL *alpha, const EREAL *a, int *lda,
   const EREAL *head_x, const EREAL *tail_x, int *incx, EREAL *beta,
   EREAL *y, int *incy, int *prec) {
  BLAS_egemv2_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *alpha,
		a, *lda, head_x, tail_x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
