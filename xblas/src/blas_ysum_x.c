#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_ysum_x(int n, const void *x, int incx,
		 void *sum, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine computes the summation:
 * 
 *     sum <- SUM_{i=0, n-1} x[i].
 * 
 * Arguments
 * =========
 *
 * n      (input) int
 *        The length of vector x.
 * 
 * x      (input) const void*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * sum    (output) void*
 * 
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: EREAL precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than EREAL, and wider than 80-bits.
 *                            We use EREAL-EREAL in our implementation.
 *
 */
{
  static const char routine_name[] = "BLAS_ysum_x";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, xi;
      EREAL *sum_i = (EREAL *) sum;
      const EREAL *x_i = (EREAL *) x;
      EREAL x_elem[2];
      EREAL tmp[2];


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -1, n, NULL);
      if (incx == 0)
	BLAS_error(routine_name, -3, incx, NULL);

      /* Immediate return. */
      if (n <= 0) {
	sum_i[0] = sum_i[1] = 0.0;
	return;
      }



      tmp[0] = tmp[1] = 0.0;

      incx *= 2;
      if (incx < 0)
	xi = -(n - 1) * incx;
      else
	xi = 0;

      for (i = 0; i < n; i++, xi += incx) {
	x_elem[0] = x_i[xi];
	x_elem[1] = x_i[xi + 1];
	tmp[0] = tmp[0] + x_elem[0];
	tmp[1] = tmp[1] + x_elem[1];
      }
      ((EREAL *) sum)[0] = tmp[0];
      ((EREAL *) sum)[1] = tmp[1];



      break;
    }

  case blas_prec_extra:
    {
      int i, xi;
      EREAL *sum_i = (EREAL *) sum;
      const EREAL *x_i = (EREAL *) x;
      EREAL x_elem[2];
      EREAL head_tmp[2], tail_tmp[2];
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -1, n, NULL);
      if (incx == 0)
	BLAS_error(routine_name, -3, incx, NULL);

      /* Immediate return. */
      if (n <= 0) {
	sum_i[0] = sum_i[1] = 0.0;
	return;
      }

      FPU_FIX_START;

      head_tmp[0] = head_tmp[1] = tail_tmp[0] = tail_tmp[1] = 0.0;

      incx *= 2;
      if (incx < 0)
	xi = -(n - 1) * incx;
      else
	xi = 0;

      for (i = 0; i < n; i++, xi += incx) {
	x_elem[0] = x_i[xi];
	x_elem[1] = x_i[xi + 1];
	{
	  EREAL head_t, tail_t;
	  EREAL head_a, tail_a;
	  head_a = head_tmp[0];
	  tail_a = tail_tmp[0];
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL. */
	    EREAL e, t1, t2;

	    /* Knuth trick. */
	    t1 = head_a + x_elem[0];
	    e = t1 - head_a;
	    t2 = ((x_elem[0] - e) + (head_a - (t1 - e))) + tail_a;

	    /* The result is t1 + t2, after normalization. */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmp[0] = head_t;
	  tail_tmp[0] = tail_t;
	  head_a = head_tmp[1];
	  tail_a = tail_tmp[1];
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL. */
	    EREAL e, t1, t2;

	    /* Knuth trick. */
	    t1 = head_a + x_elem[1];
	    e = t1 - head_a;
	    t2 = ((x_elem[1] - e) + (head_a - (t1 - e))) + tail_a;

	    /* The result is t1 + t2, after normalization. */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmp[1] = head_t;
	  tail_tmp[1] = tail_t;
	}
      }
      ((EREAL *) sum)[0] = head_tmp[0];
      ((EREAL *) sum)[1] = head_tmp[1];

      FPU_FIX_STOP;
    }
    break;
  }

}
