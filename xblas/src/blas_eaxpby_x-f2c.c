
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_eaxpby_x(int n, EREAL alpha, const EREAL *x, int incx,
		   EREAL beta, EREAL *y,
		   int incy, enum blas_prec_type prec);


extern void FC_FUNC_(blas_eaxpby_x, BLAS_EAXPBY_X)
 
  (int *n, EREAL *alpha, const EREAL *x, int *incx, EREAL *beta, EREAL *y,
   int *incy, int *prec) {
  BLAS_eaxpby_x(*n, *alpha, x, *incx, *beta, y, *incy,
		(enum blas_prec_type) *prec);
}
