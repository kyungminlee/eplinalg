
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ewaxpby_x(int n, EREAL alpha, const EREAL *x, int incx,
		    EREAL beta, const EREAL *y, int incy, EREAL *w,
		    int incw, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ewaxpby_x, BLAS_EWAXPBY_X)
 
  (int *n, EREAL *alpha, const EREAL *x, int *incx, EREAL *beta,
   const EREAL *y, int *incy, EREAL *w, int *incw, int *prec) {
  BLAS_ewaxpby_x(*n, *alpha, x, *incx, *beta, y, *incy, w, *incw,
		 (enum blas_prec_type) *prec);
}
