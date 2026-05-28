
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_egbmv2_x(enum blas_order_type order, enum blas_trans_type trans,
		   int m, int n, int kl, int ku, EREAL alpha,
		   const EREAL *a, int lda, const EREAL *head_x,
		   const EREAL *tail_x, int incx, EREAL beta,
		   EREAL *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_egbmv2_x, BLAS_EGBMV2_X)
 
  (int *trans, int *m, int *n, int *kl, int *ku, EREAL *alpha,
   const EREAL *a, int *lda, const EREAL *head_x, const EREAL *tail_x,
   int *incx, EREAL *beta, EREAL *y, int *incy, int *prec) {
  BLAS_egbmv2_x(blas_colmajor, (enum blas_trans_type) *trans, *m, *n, *kl,
		*ku, *alpha, a, *lda, head_x, tail_x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
