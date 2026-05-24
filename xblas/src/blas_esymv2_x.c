#include <blas_extended.h>
#include <blas_extended_private.h>
#include <blas_fpu.h>
void BLAS_esymv2_x(enum blas_order_type order, enum blas_uplo_type uplo,
		   int n, EREAL alpha, const EREAL *a, int lda,
		   const EREAL *x_head, const EREAL *x_tail, int incx,
		   EREAL beta, EREAL *y, int incy, enum blas_prec_type prec)

/* 
 * Purpose
 * =======
 *
 * This routines computes the matrix product:
 *
 *     y  <-  alpha * A * (x_head + x_tail) + beta * y
 * 
 * where A is a symmetric matrix.
 *
 * Arguments
 * =========
 *
 * order   (input) enum blas_order_type
 *         Storage format of input symmetric matrix A.
 * 
 * uplo    (input) enum blas_uplo_type
 *         Determines which half of matrix A (upper or lower triangle)
 *         is accessed.
 *
 * n       (input) int
 *         Dimension of A and size of vectors x, y.
 *
 * alpha   (input) EREAL
 * 
 * a       (input) EREAL*
 *         Matrix A.
 *
 * lda     (input) int
 *         Leading dimension of matrix A.
 *
 * x_head  (input) EREAL*
 *         Vector x_head
 *
 * x_tail  (input) EREAL*
 *         Vector x_tail
 *   
 * incx    (input) int
 *         Stride for vector x.
 *
 * beta    (input) EREAL
 * 
 * y       (input) EREAL*
 *         Vector y.
 *
 * incy    (input) int
 *         Stride for vector y.
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
  /* Routine name */
  const char routine_name[] = "BLAS_esymv2_x";
  switch (prec) {

  case blas_prec_single:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const EREAL *a_i = a;
      const EREAL *x_head_i = x_head;
      const EREAL *x_tail_i = x_tail;
      EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL a_elem;
      EREAL x_elem;
      EREAL y_elem;
      EREAL prod1;
      EREAL prod2;
      EREAL sum1;
      EREAL sum2;
      EREAL tmp1;
      EREAL tmp2;
      EREAL tmp3;



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	sum1 = 0.0;
	sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	sum1 = sum1 + sum2;
	tmp1 = sum1 * alpha_i;
	y_elem = y_i[yi];
	tmp2 = y_elem * beta_i;
	tmp3 = tmp1 + tmp2;
	y_i[yi] = tmp3;
      }



      break;
    }

  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const EREAL *a_i = a;
      const EREAL *x_head_i = x_head;
      const EREAL *x_tail_i = x_tail;
      EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL a_elem;
      EREAL x_elem;
      EREAL y_elem;
      EREAL prod1;
      EREAL prod2;
      EREAL sum1;
      EREAL sum2;
      EREAL tmp1;
      EREAL tmp2;
      EREAL tmp3;



      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);



      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	sum1 = 0.0;
	sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  prod1 = a_elem * x_elem;
	  sum1 = sum1 + prod1;
	  x_elem = x_tail_i[xi];
	  prod2 = a_elem * x_elem;
	  sum2 = sum2 + prod2;
	}
	sum1 = sum1 + sum2;
	tmp1 = sum1 * alpha_i;
	y_elem = y_i[yi];
	tmp2 = y_elem * beta_i;
	tmp3 = tmp1 + tmp2;
	y_i[yi] = tmp3;
      }



      break;
    }

  case blas_prec_extra:{

      int i, j;
      int xi, yi, xi0, yi0;
      int aij, ai;
      int incai;
      int incaij, incaij2;

      const EREAL *a_i = a;
      const EREAL *x_head_i = x_head;
      const EREAL *x_tail_i = x_tail;
      EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL a_elem;
      EREAL x_elem;
      EREAL y_elem;
      EREAL head_prod1, tail_prod1;
      EREAL head_prod2, tail_prod2;
      EREAL head_sum1, tail_sum1;
      EREAL head_sum2, tail_sum2;
      EREAL head_tmp1, tail_tmp1;
      EREAL head_tmp2, tail_tmp2;
      EREAL head_tmp3, tail_tmp3;

      FPU_FIX_DECL;

      /* Test for no-op */
      if (n <= 0) {
	return;
      }
      if (alpha_i == 0.0 && beta_i == 1.0) {
	return;
      }

      /* Check for error conditions. */
      if (n < 0) {
	BLAS_error(routine_name, -3, n, NULL);
      }
      if (lda < n) {
	BLAS_error(routine_name, -6, n, NULL);
      }
      if (incx == 0) {
	BLAS_error(routine_name, -9, incx, NULL);
      }
      if (incy == 0) {
	BLAS_error(routine_name, -12, incy, NULL);
      }

      if ((order == blas_colmajor && uplo == blas_upper) ||
	  (order == blas_rowmajor && uplo == blas_lower)) {
	incai = lda;
	incaij = 1;
	incaij2 = lda;
      } else {
	incai = 1;
	incaij = lda;
	incaij2 = 1;
      }






      xi0 = (incx > 0) ? 0 : ((-n + 1) * incx);
      yi0 = (incy > 0) ? 0 : ((-n + 1) * incy);

      FPU_FIX_START;

      /* The most general form,   y <--- alpha * A * (x_head + x_tail) + beta * y   */
      for (i = 0, yi = yi0, ai = 0; i < n; i++, yi += incy, ai += incai) {
	head_sum1 = tail_sum1 = 0.0;
	head_sum2 = tail_sum2 = 0.0;

	for (j = 0, aij = ai, xi = xi0; j < i; j++, aij += incaij, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = a_elem * split;
	    a1 = con - a_elem;
	    a1 = con - a1;
	    a2 = a_elem - a1;
	    con = x_elem * split;
	    b1 = con - x_elem;
	    b1 = con - b1;
	    b2 = x_elem - b1;

	    head_prod1 = a_elem * x_elem;
	    tail_prod1 =
	      (((a1 * b1 - head_prod1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum1 + head_prod1;
	    bv = s1 - head_sum1;
	    s2 = ((head_prod1 - bv) + (head_sum1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum1 + tail_prod1;
	    bv = t1 - tail_sum1;
	    t2 = ((tail_prod1 - bv) + (tail_sum1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum1 = t1 + t2;
	    tail_sum1 = t2 - (head_sum1 - t1);
	  }
	  x_elem = x_tail_i[xi];
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = a_elem * split;
	    a1 = con - a_elem;
	    a1 = con - a1;
	    a2 = a_elem - a1;
	    con = x_elem * split;
	    b1 = con - x_elem;
	    b1 = con - b1;
	    b2 = x_elem - b1;

	    head_prod2 = a_elem * x_elem;
	    tail_prod2 =
	      (((a1 * b1 - head_prod2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum2 + head_prod2;
	    bv = s1 - head_sum2;
	    s2 = ((head_prod2 - bv) + (head_sum2 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum2 + tail_prod2;
	    bv = t1 - tail_sum2;
	    t2 = ((tail_prod2 - bv) + (tail_sum2 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum2 = t1 + t2;
	    tail_sum2 = t2 - (head_sum2 - t1);
	  }
	}
	for (; j < n; j++, aij += incaij2, xi += incx) {
	  a_elem = a_i[aij];
	  x_elem = x_head_i[xi];
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = a_elem * split;
	    a1 = con - a_elem;
	    a1 = con - a1;
	    a2 = a_elem - a1;
	    con = x_elem * split;
	    b1 = con - x_elem;
	    b1 = con - b1;
	    b2 = x_elem - b1;

	    head_prod1 = a_elem * x_elem;
	    tail_prod1 =
	      (((a1 * b1 - head_prod1) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum1 + head_prod1;
	    bv = s1 - head_sum1;
	    s2 = ((head_prod1 - bv) + (head_sum1 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum1 + tail_prod1;
	    bv = t1 - tail_sum1;
	    t2 = ((tail_prod1 - bv) + (tail_sum1 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum1 = t1 + t2;
	    tail_sum1 = t2 - (head_sum1 - t1);
	  }
	  x_elem = x_tail_i[xi];
	  {
	    /* Compute double_double = EREAL * EREAL. */
	    EREAL a1, a2, b1, b2, con;

	    con = a_elem * split;
	    a1 = con - a_elem;
	    a1 = con - a1;
	    a2 = a_elem - a1;
	    con = x_elem * split;
	    b1 = con - x_elem;
	    b1 = con - b1;
	    b2 = x_elem - b1;

	    head_prod2 = a_elem * x_elem;
	    tail_prod2 =
	      (((a1 * b1 - head_prod2) + a1 * b2) + a2 * b1) + a2 * b2;
	  }
	  {
	    /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	    EREAL bv;
	    EREAL s1, s2, t1, t2;

	    /* Add two hi words. */
	    s1 = head_sum2 + head_prod2;
	    bv = s1 - head_sum2;
	    s2 = ((head_prod2 - bv) + (head_sum2 - (s1 - bv)));

	    /* Add two lo words. */
	    t1 = tail_sum2 + tail_prod2;
	    bv = t1 - tail_sum2;
	    t2 = ((tail_prod2 - bv) + (tail_sum2 - (t1 - bv)));

	    s2 += t1;

	    /* Renormalize (s1, s2)  to  (t1, s2) */
	    t1 = s1 + s2;
	    s2 = s2 - (t1 - s1);

	    t2 += s2;

	    /* Renormalize (t1, t2)  */
	    head_sum2 = t1 + t2;
	    tail_sum2 = t2 - (head_sum2 - t1);
	  }
	}
	{
	  /* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
	  EREAL bv;
	  EREAL s1, s2, t1, t2;

	  /* Add two hi words. */
	  s1 = head_sum1 + head_sum2;
	  bv = s1 - head_sum1;
	  s2 = ((head_sum2 - bv) + (head_sum1 - (s1 - bv)));

	  /* Add two lo words. */
	  t1 = tail_sum1 + tail_sum2;
	  bv = t1 - tail_sum1;
	  t2 = ((tail_sum2 - bv) + (tail_sum1 - (t1 - bv)));

	  s2 += t1;

	  /* Renormalize (s1, s2)  to  (t1, s2) */
	  t1 = s1 + s2;
	  s2 = s2 - (t1 - s1);

	  t2 += s2;

	  /* Renormalize (t1, t2)  */
	  head_sum1 = t1 + t2;
	  tail_sum1 = t2 - (head_sum1 - t1);
	}
	{
	  /* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
	  EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

	  con = head_sum1 * split;
	  a11 = con - head_sum1;
	  a11 = con - a11;
	  a21 = head_sum1 - a11;
	  con = alpha_i * split;
	  b1 = con - alpha_i;
	  b1 = con - b1;
	  b2 = alpha_i - b1;

	  c11 = head_sum1 * alpha_i;
	  c21 = (((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

	  c2 = tail_sum1 * alpha_i;
	  t1 = c11 + c2;
	  t2 = (c2 - (t1 - c11)) + c21;

	  head_tmp1 = t1 + t2;
	  tail_tmp1 = t2 - (head_tmp1 - t1);
	}
	y_elem = y_i[yi];
	{
	  /* Compute double_double = EREAL * EREAL. */
	  EREAL a1, a2, b1, b2, con;

	  con = y_elem * split;
	  a1 = con - y_elem;
	  a1 = con - a1;
	  a2 = y_elem - a1;
	  con = beta_i * split;
	  b1 = con - beta_i;
	  b1 = con - b1;
	  b2 = beta_i - b1;

	  head_tmp2 = y_elem * beta_i;
	  tail_tmp2 = (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	}
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
	  head_tmp3 = t1 + t2;
	  tail_tmp3 = t2 - (head_tmp3 - t1);
	}
	y_i[yi] = head_tmp3;
      }

      FPU_FIX_STOP;

      break;
    }
  }
}				/* end BLAS_esymv2_x */
