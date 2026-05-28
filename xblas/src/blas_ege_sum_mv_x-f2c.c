
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ege_sum_mv_x(enum blas_order_type order, int m, int n,
		       EREAL alpha, const EREAL *a, int lda,
		       const EREAL *x, int incx,
		       EREAL beta, const EREAL *b, int ldb,
		       EREAL *y, int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ege_sum_mv_x, BLAS_EGE_SUM_MV_X)
 
  (int *m, int *n, EREAL *alpha, const EREAL *a, int *lda, const EREAL *x,
   int *incx, EREAL *beta, const EREAL *b, int *ldb, EREAL *y, int *incy,
   int *prec) {
  BLAS_ege_sum_mv_x(blas_colmajor, *m, *n, *alpha, a, *lda, x, *incx, *beta,
		    b, *ldb, y, *incy, (enum blas_prec_type) *prec);
}
