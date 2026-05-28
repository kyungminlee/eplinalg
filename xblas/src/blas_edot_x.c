#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_edot_x(enum blas_conj_type conj, int n, EREAL alpha,
		 const EREAL *x, int incx, EREAL beta,
		 const EREAL *y, int incy,
		 EREAL *r, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 * 
 * This routine computes the inner product:
 * 
 *     r <- beta * r + alpha * SUM_{i=0, n-1} x[i] * y[i].
 * 
 * Arguments
 * =========
 *  
 * conj   (input) enum blas_conj_type
 *        When x and y are complex vectors, specifies whether vector
 *        components x[i] are used unconjugated or conjugated. 
 * 
 * n      (input) int
 *        The length of vectors x and y.
 * 
 * alpha  (input) EREAL
 * 
 * x      (input) const EREAL*
 *        Array of length n.
 * 
 * incx   (input) int
 *        The stride used to access components x[i].
 *
 * beta   (input) EREAL
 *
 * y      (input) const EREAL*
 *        Array of length n.
 *      
 * incy   (input) int
 *        The stride used to access components y[i].
 *
 * r      (input/output) EREAL*
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
  static const char routine_name[] = "BLAS_edot_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      EREAL *r_i = r;
      const EREAL *x_i = x;
      const EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL x_ii;
      EREAL y_ii;
      EREAL r_v;
      EREAL prod;
      EREAL sum;
      EREAL tmp1;
      EREAL tmp2;


      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if ((beta_i == 1.0) && (n == 0 || (alpha_i == 0.0)))
	return;



      r_v = r_i[0];
      sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	prod = x_ii * y_ii;	/* prod = x[i]*y[i] */
	sum = sum + prod;	/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      tmp1 = sum * alpha_i;	/* tmp1 = sum*alpha */
      tmp2 = r_v * beta_i;	/* tmp2 = r*beta */
      tmp1 = tmp1 + tmp2;	/* tmp1 = tmp1+tmp2 */
      *r = tmp1;		/* r = tmp1 */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      EREAL *r_i = r;
      const EREAL *x_i = x;
      const EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL x_ii;
      EREAL y_ii;
      EREAL r_v;
      EREAL head_prod, tail_prod;
      EREAL head_sum, tail_sum;
      EREAL head_tmp1, tail_tmp1;
      EREAL head_tmp2, tail_tmp2;
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (n < 0)
	BLAS_error(routine_name, -2, n, NULL);
      else if (incx == 0)
	BLAS_error(routine_name, -5, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -8, incy, NULL);

      /* Immediate return. */
      if ((beta_i == 1.0) && (n == 0 || (alpha_i == 0.0)))
	return;

      FPU_FIX_START;

      r_v = r_i[0];
      head_sum = tail_sum = 0.0;


      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii = x_i[ix];
	y_ii = y_i[iy];

	{
	  /* Compute double_double = EREAL * EREAL. */
	  EREAL a1, a2, b1, b2, con;

	  con = x_ii * split;
	  a1 = con - x_ii;
	  a1 = con - a1;
	  a2 = x_ii - a1;
	  con = y_ii * split;
	  b1 = con - y_ii;
	  b1 = con - b1;
	  b2 = y_ii - b1;

	  head_prod = x_ii * y_ii;
	  tail_prod = (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	}			/* prod = x[i]*y[i] */
	{
	  /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	  EREAL bv;
	  EREAL s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_sum + head_prod;
	  bv = s1 - head_sum;
	  s2 = ((head_prod - bv) + (head_sum - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_sum + tail_prod;
	  bv = t1 - tail_sum;
	  t2 = ((tail_prod - bv) + (tail_sum - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_sum = t1 + t2;
	  tail_sum = t2 - (head_sum - t1);
	}			/* sum = sum+prod */
	ix += incx;
	iy += incy;
      }				/* endfor */


      {
	/* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
	EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	con = head_sum * split;
	a11 = con - head_sum;
	a11 = con - a11;
	a21 = head_sum - a11;
	con = alpha_i * split;
	b1 = con - alpha_i;
	b1 = con - b1;
	b2 = alpha_i - b1;

	c11 = head_sum * alpha_i;
	c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	c2 = tail_sum * alpha_i;
	t1 = c11 + c2;
	t2 = (c2 - (t1 - c11)) + c21;

	head_tmp1 = t1 + t2;
	tail_tmp1 = t2 - (head_tmp1 - t1);
      }				/* tmp1 = sum*alpha */
      {
	/* Compute double_double = EREAL * EREAL. */
	EREAL a1, a2, b1, b2, con;

	con = r_v * split;
	a1 = con - r_v;
	a1 = con - a1;
	a2 = r_v - a1;
	con = beta_i * split;
	b1 = con - beta_i;
	b1 = con - b1;
	b2 = beta_i - b1;

	head_tmp2 = r_v * beta_i;
	tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
      }				/* tmp2 = r*beta */
      {
	/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	EREAL bv;
	EREAL s1, s2, t1, t2;

	/* Add two hi words. */
	s1 = head_tmp1 + head_tmp2;
	bv = s1 - head_tmp1;
	s2 = ((head_tmp2 - bv) + (head_tmp1 - (s1 - bv)));

	/* Add two lo words. */
	t1 = tail_tmp1 + tail_tmp2;
	bv = t1 - tail_tmp1;
	t2 = ((tail_tmp2 - bv) + (tail_tmp1 - (t1 - bv)));

	s2 += t1;

	/* Renormalize (s1, s2)  to  (t1, s2) */
	t1 = s1 + s2;
	s2 = s2 - (t1 - s1);

	t2 += s2;

	/* Renormalize (t1, t2)  */
	head_tmp1 = t1 + t2;
	tail_tmp1 = t2 - (head_tmp1 - t1);
      }				/* tmp1 = tmp1+tmp2 */
      *r = head_tmp1;		/* r = tmp1 */

      FPU_FIX_STOP;
    }
    break;
  }
}
