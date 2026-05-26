
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_ysum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec);


extern void FC_FUNC_(blas_ysum_x, BLAS_YSUM_X)
  (int *n, const void *x, int *incx, void *sum, int *prec) {
  BLAS_ysum_x(*n, x, *incx, sum, (enum blas_prec_type) *prec);
}
