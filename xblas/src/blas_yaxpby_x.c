#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_yaxpby_x(int n, const void *alpha, const void *x, int incx,
		   const void *beta, void *y,
		   int incy, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * This routine computes:
 *
 *      y <- alpha * x + beta * y.
 *
 * Arguments
 * =========
 * 
 * n         (input) int
 *           The length of vectors x and y.
 * 
 * alpha     (input) const void*
 *
 * x         (input) const void*
 *           Array of length n.
 *
 * incx      (input) int
 *           The stride used to access components x[i].
 * 
 * beta      (input) const void*
 *
 * y         (input) void*
 *           Array of length n.
 * 
 * incy      (input) int
 *           The stride used to access components y[i].
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
  static const char routine_name[] = "BLAS_yaxpby_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:
    {
      int i, ix = 0, iy = 0;
      const EREAL *x_i = (EREAL *) x;
      EREAL *y_i = (EREAL *) y;
      EREAL *alpha_i = (EREAL *) alpha;
      EREAL *beta_i = (EREAL *) beta;
      EREAL x_ii[2];
      EREAL y_ii[2];
      EREAL tmpx[2];
      EREAL tmpy[2];


      /* Test the input parameters. */
      if (incx == 0)
	BLAS_error(routine_name, -4, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -7, incy, NULL);

      /* Immediate return */
      if (n <= 0
	  || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	      && (beta_i[0] == 1.0 && beta_i[1] == 0.0)))
	return;



      incx *= 2;
      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii[0] = x_i[ix];
	x_ii[1] = x_i[ix + 1];
	y_ii[0] = y_i[iy];
	y_ii[1] = y_i[iy + 1];
	{
	  tmpx[0] =
	    (EREAL) alpha_i[0] * x_ii[0] - (EREAL) alpha_i[1] * x_ii[1];
	  tmpx[1] =
	    (EREAL) alpha_i[0] * x_ii[1] + (EREAL) alpha_i[1] * x_ii[0];
	}			/* tmpx  = alpha * x[ix] */
	{
	  tmpy[0] =
	    (EREAL) beta_i[0] * y_ii[0] - (EREAL) beta_i[1] * y_ii[1];
	  tmpy[1] =
	    (EREAL) beta_i[0] * y_ii[1] + (EREAL) beta_i[1] * y_ii[0];
	}			/* tmpy = beta * y[iy] */
	tmpy[0] = tmpy[0] + tmpx[0];
	tmpy[1] = tmpy[1] + tmpx[1];
	y_i[iy] = tmpy[0];
	y_i[iy + 1] = tmpy[1];
	ix += incx;
	iy += incy;
      }				/* endfor */


    }
    break;
  case blas_prec_extra:
    {
      int i, ix = 0, iy = 0;
      const EREAL *x_i = (EREAL *) x;
      EREAL *y_i = (EREAL *) y;
      EREAL *alpha_i = (EREAL *) alpha;
      EREAL *beta_i = (EREAL *) beta;
      EREAL x_ii[2];
      EREAL y_ii[2];
      EREAL head_tmpx[2], tail_tmpx[2];
      EREAL head_tmpy[2], tail_tmpy[2];
      FPU_FIX_DECL;

      /* Test the input parameters. */
      if (incx == 0)
	BLAS_error(routine_name, -4, incx, NULL);
      else if (incy == 0)
	BLAS_error(routine_name, -7, incy, NULL);

      /* Immediate return */
      if (n <= 0
	  || (alpha_i[0] == 0.0 && alpha_i[1] == 0.0
	      && (beta_i[0] == 1.0 && beta_i[1] == 0.0)))
	return;

      FPU_FIX_START;

      incx *= 2;
      incy *= 2;
      if (incx < 0)
	ix = (-n + 1) * incx;
      if (incy < 0)
	iy = (-n + 1) * incy;

      for (i = 0; i < n; ++i) {
	x_ii[0] = x_i[ix];
	x_ii[1] = x_i[ix + 1];
	y_ii[0] = y_i[iy];
	y_ii[1] = y_i[iy + 1];
	{
	  /* Compute complex-extra = complex-EREAL * complex-EREAL. */
	  EREAL head_t1, tail_t1;
	  EREAL head_t2, tail_t2;
	  /* Real part */
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = alpha_i[0] * split;
	    a1 = con - alpha_i[0];
	    a1 = con - a1;
	    a2 = alpha_i[0] - a1;
	    con = x_ii[0] * split;
	    b1 = con - x_ii[0];
	    b1 = con - b1;
	    b2 = x_ii[0] - b1;

	    head_t1 = alpha_i[0] * x_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = alpha_i[1] * split;
	    a1 = con - alpha_i[1];
	    a1 = con - a1;
	    a2 = alpha_i[1] - a1;
	    con = x_ii[1] * split;
	    b1 = con - x_ii[1];
	    b1 = con - b1;
	    b2 = x_ii[1] - b1;

	    head_t2 = alpha_i[1] * x_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  head_t2 = -head_t2;
	  tail_t2 = -tail_t2;
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpx[0] = head_t1;
	  tail_tmpx[0] = tail_t1;
	  /* Imaginary part */
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = alpha_i[1] * split;
	    a1 = con - alpha_i[1];
	    a1 = con - a1;
	    a2 = alpha_i[1] - a1;
	    con = x_ii[0] * split;
	    b1 = con - x_ii[0];
	    b1 = con - b1;
	    b2 = x_ii[0] - b1;

	    head_t1 = alpha_i[1] * x_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = alpha_i[0] * split;
	    a1 = con - alpha_i[0];
	    a1 = con - a1;
	    a2 = alpha_i[0] - a1;
	    con = x_ii[1] * split;
	    b1 = con - x_ii[1];
	    b1 = con - b1;
	    b2 = x_ii[1] - b1;

	    head_t2 = alpha_i[0] * x_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpx[1] = head_t1;
	  tail_tmpx[1] = tail_t1;
	}			/* tmpx  = alpha * x[ix] */
	{
	  /* Compute complex-extra = complex-EREAL * complex-EREAL. */
	  EREAL head_t1, tail_t1;
	  EREAL head_t2, tail_t2;
	  /* Real part */
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_ii[0] * split;
	    b1 = con - y_ii[0];
	    b1 = con - b1;
	    b2 = y_ii[0] - b1;

	    head_t1 = beta_i[0] * y_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_ii[1] * split;
	    b1 = con - y_ii[1];
	    b1 = con - b1;
	    b2 = y_ii[1] - b1;

	    head_t2 = beta_i[1] * y_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  head_t2 = -head_t2;
	  tail_t2 = -tail_t2;
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpy[0] = head_t1;
	  tail_tmpy[0] = tail_t1;
	  /* Imaginary part */
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = beta_i[1] * split;
	    a1 = con - beta_i[1];
	    a1 = con - a1;
	    a2 = beta_i[1] - a1;
	    con = y_ii[0] * split;
	    b1 = con - y_ii[0];
	    b1 = con - b1;
	    b2 = y_ii[0] - b1;

	    head_t1 = beta_i[1] * y_ii[0];
	    tail_t1 = (((a1 * b1 - head_t1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = beta_i[0] * split;
	    a1 = con - beta_i[0];
	    a1 = con - a1;
	    a2 = beta_i[0] - a1;
	    con = y_ii[1] * split;
	    b1 = con - y_ii[1];
	    b1 = con - b1;
	    b2 = y_ii[1] - b1;

	    head_t2 = beta_i[0] * y_ii[1];
	    tail_t2 = (((a1 * b1 - head_t2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_t1 + head_t2;
	    bv = s1 - head_t1;
	    s2 = ((head_t2 - bv) + (head_t1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_t1 + tail_t2;
	    bv = t1 - tail_t1;
	    t2 = ((tail_t2 - bv) + (tail_t1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t1 = t1 + t2;
	    tail_t1 = t2 - (head_t1 - t1);
	  }
	  head_tmpy[1] = head_t1;
	  tail_tmpy[1] = tail_t1;
	}			/* tmpy = beta * y[iy] */
	{
	  EREAL head_t, tail_t;
	  EREAL head_a, tail_a;
	  EREAL head_b, tail_b;
	  /* Real part */
	  head_a = head_tmpy[0];
	  tail_a = tail_tmpy[0];
	  head_b = head_tmpx[0];
	  tail_b = tail_tmpx[0];
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_a + head_b;
	    bv = s1 - head_a;
	    s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_a + tail_b;
	    bv = t1 - tail_a;
	    t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmpy[0] = head_t;
	  tail_tmpy[0] = tail_t;
	  /* Imaginary part */
	  head_a = head_tmpy[1];
	  tail_a = tail_tmpy[1];
	  head_b = head_tmpx[1];
	  tail_b = tail_tmpx[1];
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_a + head_b;
	    bv = s1 - head_a;
	    s2 = ((head_b - bv) + (head_a - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_a + tail_b;
	    bv = t1 - tail_a;
	    t2 = ((tail_b - bv) + (tail_a - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_t = t1 + t2;
	    tail_t = t2 - (head_t - t1);
	  }
	  head_tmpy[1] = head_t;
	  tail_tmpy[1] = tail_t;
	}
	y_i[iy] = head_tmpy[0];
	y_i[iy + 1] = head_tmpy[1];
	ix += incx;
	iy += incy;
      }				/* endfor */

      FPU_FIX_STOP;
    }
    break;
  }
}				/* end BLAS_yaxpby_x */
