
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_esum_x(int n, const EREAL *x, int incx,
		 EREAL *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_esum_x, BLAS_ESUM_X)
  (int *n, const EREAL *x, int *incx, EREAL *sum, int *prec) {
  BLAS_esum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
