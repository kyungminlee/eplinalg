
#include "f2c-bridge.h"
#include "blas_enum.h"
void BLAS_edot_x(enum blas_conj_type conj, int n, EREAL alpha,
		 const EREAL *x, int incx, EREAL beta,
		 const EREAL *y, int incy,
		 EREAL *r, enum blas_prec_type prec);


extern void FC_FUNC_(blas_edot_x, BLAS_EDOT_X)
 
  (int *conj, int *n, EREAL *alpha, const EREAL *x, int *incx, EREAL *beta,
   const EREAL *y, int *incy, EREAL *r, int *prec) {
  BLAS_edot_x((enum blas_conj_type) *conj, *n, *alpha, x, *incx, *beta, y,
	      *incy, r, (enum blas_prec_type) *prec);
}
