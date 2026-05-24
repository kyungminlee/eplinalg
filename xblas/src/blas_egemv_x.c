#include "blas_extended.h"
#include "blas_extended_private.h"
void BLAS_egemv_x(enum blas_order_type order, enum blas_trans_type trans,
		  int m, int n, EREAL alpha, const EREAL *a, int lda,
		  const EREAL *x, int incx, EREAL beta, EREAL *y,
		  int incy, enum blas_prec_type prec)

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * A * x + beta * y, where A is a general matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of AP; row or column major
 *
 * trans        (input) blas_trans_type
 *              Transpose of AB; no trans, 
 *              trans, or conjugate trans
 *
 * m            (input) int
 *              Dimension of AB
 *
 * n            (input) int
 *              Dimension of AB and the length of vector x
 *
 * alpha        (input) EREAL
 *              
 * A            (input) const EREAL*
 *
 * lda          (input) int 
 *              Leading dimension of A
 *
 * x            (input) const EREAL*
 * 
 * incx         (input) int
 *              The stride for vector x.
 *
 * beta         (input) EREAL
 *
 * y            (input/output) EREAL*
 *
 * incy         (input) int
 *              The stride for vector y.
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
  static const char routine_name[] = "BLAS_egemv_x";
  switch (prec) {
  case blas_prec_single:
  case blas_prec_double:
  case blas_prec_indigenous:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const EREAL *a_i = a;
      const EREAL *x_i = x;
      EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL a_elem;
      EREAL x_elem;
      EREAL y_elem;
      EREAL prod;
      EREAL sum;
      EREAL tmp1;
      EREAL tmp2;


      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -9, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -12, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }
      if ((order == blas_colmajor && lda < m) ||
	  (order == blas_rowmajor && lda < n))
	BLAS_error(routine_name, -7, lda, NULL);








      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem = y_i[iy];
	    tmp1 = y_elem * beta_i;
	    y_i[iy] = tmp1;
	    iy += incy;
	  }
	}
      } else {

	/* if beta = 0, we can save m multiplies: y = alpha*A*x */
	if (beta_i == 0.0) {
	  /* save m more multiplies if alpha = 1 */
	  if (alpha_i == 1.0) {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		aij += incaij;
		jx += incx;
	      }
	      y_i[iy] = sum;
	      ai += incai;
	      iy += incy;
	    }
	  } else {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
		prod = a_elem * x_elem;
		sum = sum + prod;
		aij += incaij;
		jx += incx;
	      }
	      tmp1 = sum * alpha_i;
	      y_i[iy] = tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {
	  /* the most general form, y = alpha*A*x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    sum = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem = a_i[aij];

	      x_elem = x_i[jx];
	      prod = a_elem * x_elem;
	      sum = sum + prod;
	      aij += incaij;
	      jx += incx;
	    }
	    tmp1 = sum * alpha_i;
	    y_elem = y_i[iy];
	    tmp2 = y_elem * beta_i;
	    tmp1 = tmp1 + tmp2;
	    y_i[iy] = tmp1;
	    ai += incai;
	    iy += incy;
	  }
	}

      }



      break;
    }
  case blas_prec_extra:{

      int i, j;
      int iy, jx, kx, ky;
      int lenx, leny;
      int ai, aij;
      int incai, incaij;

      const EREAL *a_i = a;
      const EREAL *x_i = x;
      EREAL *y_i = y;
      EREAL alpha_i = alpha;
      EREAL beta_i = beta;
      EREAL a_elem;
      EREAL x_elem;
      EREAL y_elem;
      EREAL head_prod, tail_prod;
      EREAL head_sum, tail_sum;
      EREAL head_tmp1, tail_tmp1;
      EREAL head_tmp2, tail_tmp2;
      FPU_FIX_DECL;

      /* all error calls */
      if (m < 0)
	BLAS_error(routine_name, -3, m, 0);
      else if (n <= 0)
	BLAS_error(routine_name, -4, n, 0);
      else if (incx == 0)
	BLAS_error(routine_name, -9, incx, 0);
      else if (incy == 0)
	BLAS_error(routine_name, -12, incy, 0);

      if ((order == blas_rowmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = lda;
	incaij = 1;
      } else if ((order == blas_rowmajor) && (trans != blas_no_trans)) {
	lenx = m;
	leny = n;
	incai = 1;
	incaij = lda;
      } else if ((order == blas_colmajor) && (trans == blas_no_trans)) {
	lenx = n;
	leny = m;
	incai = 1;
	incaij = lda;
      } else {			/* colmajor and blas_trans */
	lenx = m;
	leny = n;
	incai = lda;
	incaij = 1;
      }
      if ((order == blas_colmajor && lda < m) ||
	  (order == blas_rowmajor && lda < n))
	BLAS_error(routine_name, -7, lda, NULL);

      FPU_FIX_START;






      if (incx > 0)
	kx = 0;
      else
	kx = (1 - lenx) * incx;
      if (incy > 0)
	ky = 0;
      else
	ky = (1 - leny) * incy;

      /* No extra-precision needed for alpha = 0 */
      if (alpha_i == 0.0) {
	if (beta_i == 0.0) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_i[iy] = 0.0;
	    iy += incy;
	  }
	} else if (!(beta_i == 0.0)) {
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    y_elem = y_i[iy];
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

	      head_tmp1 = y_elem * beta_i;
	      tail_tmp1 =
		(((a1 * b1 - head_tmp1) + a1 * b2) + a2 * b1) + a2 * b2;
	    }
	    y_i[iy] = head_tmp1;
	    iy += incy;
	  }
	}
      } else {

	/* if beta = 0, we can save m multiplies: y = alpha*A*x */
	if (beta_i == 0.0) {
	  /* save m more multiplies if alpha = 1 */
	  if (alpha_i == 1.0) {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
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

		  head_prod = a_elem * x_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
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
		}
		aij += incaij;
		jx += incx;
	      }
	      y_i[iy] = head_sum;
	      ai += incai;
	      iy += incy;
	    }
	  } else {
	    ai = 0;
	    iy = ky;
	    for (i = 0; i < leny; i++) {
	      head_sum = tail_sum = 0.0;
	      aij = ai;
	      jx = kx;
	      for (j = 0; j < lenx; j++) {
		a_elem = a_i[aij];

		x_elem = x_i[jx];
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

		  head_prod = a_elem * x_elem;
		  tail_prod =
		    (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
		}
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
		}
		aij += incaij;
		jx += incx;
	      }
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
	      }
	      y_i[iy] = head_tmp1;
	      ai += incai;
	      iy += incy;
	    }
	  }
	} else {
	  /* the most general form, y = alpha*A*x + beta*y */
	  ai = 0;
	  iy = ky;
	  for (i = 0; i < leny; i++) {
	    head_sum = tail_sum = 0.0;;
	    aij = ai;
	    jx = kx;
	    for (j = 0; j < lenx; j++) {
	      a_elem = a_i[aij];

	      x_elem = x_i[jx];
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

		head_prod = a_elem * x_elem;
		tail_prod =
		  (((a1 * b1 - head_prod) + a1 * b2) + a2 * b1) + a2 * b2;
	      }
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
	      }
	      aij += incaij;
	      jx += incx;
	    }
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
	    }
	    y_elem = y_i[iy];
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
	      tail_tmp2 =
		(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
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
	      head_tmp1 = t1 + t2;
	      tail_tmp1 = t2 - (head_tmp1 - t1);
	    }
	    y_i[iy] = head_tmp1;
	    ai += incai;
	    iy += incy;
	  }
	}

      }

      FPU_FIX_STOP;
    }
    break;
  }
}
