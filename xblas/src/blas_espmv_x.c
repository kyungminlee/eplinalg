#include "blas_extended.h"
#include "blas_extended_private.h"

/*
 * Purpose
 * =======
 *
 * Computes y = alpha * ap * x + beta * y, where ap is a symmetric
 * packed matrix.
 *
 * Arguments
 * =========
 *
 * order        (input) blas_order_type
 *              Order of ap; row or column major
 *
 * uplo         (input) blas_uplo_type
 *              Whether ap is upper or lower
 *
 * n            (input) int
 *              Dimension of ap and the length of vector x
 *
 * alpha        (input) EREAL
 *              
 * ap           (input) EREAL*
 *
 * x            (input) EREAL*
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
 * prec   (input) enum blas_prec_type
 *        Specifies the internal precision to be used.
 *        = blas_prec_single: single precision.
 *        = blas_prec_double: EREAL precision.
 *        = blas_prec_extra : anything at least 1.5 times as accurate
 *                            than EREAL, and wider than 80-bits.
 *                            We use EREAL-EREAL in our implementation.
 *
 *
 */
void BLAS_espmv_x(enum blas_order_type order, enum blas_uplo_type uplo,
		  int n, EREAL alpha, const EREAL *ap,
		  const EREAL *x, int incx, EREAL beta,
		  EREAL *y, int incy, enum blas_prec_type prec)
{
  static const char routine_name[] = "BLAS_espmv_x";

  switch (prec) {
  case blas_prec_single:
  case blas_prec_indigenous:
  case blas_prec_double:{
      {
	int matrix_row, step, ap_index, ap_start, x_index, x_start;
	int y_start, y_index, incap;
	EREAL alpha_i = alpha;
	EREAL beta_i = beta;

	const EREAL *ap_i = ap;
	const EREAL *x_i = x;
	EREAL *y_i = y;
	EREAL rowsum;
	EREAL rowtmp;
	EREAL matval;
	EREAL vecval;
	EREAL resval;
	EREAL tmp1;
	EREAL tmp2;


	incap = 1;




	if (incx < 0)
	  x_start = (-n + 1) * incx;
	else
	  x_start = 0;
	if (incy < 0)
	  y_start = (-n + 1) * incy;
	else
	  y_start = 0;

	if (n < 1) {
	  return;
	}
	if (alpha_i == 0.0 && beta_i == 1.0) {
	  return;
	}

	/* Check for error conditions. */
	if (order != blas_colmajor && order != blas_rowmajor) {
	  BLAS_error(routine_name, -1, order, NULL);
	}
	if (uplo != blas_upper && uplo != blas_lower) {
	  BLAS_error(routine_name, -2, uplo, NULL);
	}
	if (incx == 0) {
	  BLAS_error(routine_name, -7, incx, NULL);
	}
	if (incy == 0) {
	  BLAS_error(routine_name, -10, incy, NULL);
	}



	if (alpha_i == 0.0) {
	  {
	    y_index = y_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      resval = y_i[y_index];

	      tmp2 = beta_i * resval;

	      y_i[y_index] = tmp2;

	      y_index += incy;
	    }
	  }
	} else {
	  if (uplo == blas_lower)
	    order = (order == blas_rowmajor) ? blas_colmajor : blas_rowmajor;
	  if (order == blas_rowmajor) {
	    if (alpha_i == 1.0) {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    tmp1 = rowsum;
		    y_i[y_index] = tmp1;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    tmp1 = rowsum;
		    tmp2 = beta_i * resval;
		    tmp2 = tmp1 + tmp2;
		    y_i[y_index] = tmp2;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      }
	    } else {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    tmp1 = rowsum * alpha_i;
		    y_i[y_index] = tmp1;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    tmp1 = rowsum * alpha_i;
		    tmp2 = beta_i * resval;
		    tmp2 = tmp1 + tmp2;
		    y_i[y_index] = tmp2;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      }
	    }
	  } else {
	    if (alpha_i == 1.0) {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    tmp1 = rowsum;
		    y_i[y_index] = tmp1;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    tmp1 = rowsum;
		    tmp2 = beta_i * resval;
		    tmp2 = tmp1 + tmp2;
		    y_i[y_index] = tmp2;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      }
	    } else {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    tmp1 = rowsum * alpha_i;
		    y_i[y_index] = tmp1;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    rowsum = 0.0;
		    rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      rowtmp = matval * vecval;
		      rowsum = rowsum + rowtmp;
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    tmp1 = rowsum * alpha_i;
		    tmp2 = beta_i * resval;
		    tmp2 = tmp1 + tmp2;
		    y_i[y_index] = tmp2;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      }
	    }
	  }			/* if order == ... */
	}			/* alpha != 0 */


      }
      break;
    }

  case blas_prec_extra:{
      {
	int matrix_row, step, ap_index, ap_start, x_index, x_start;
	int y_start, y_index, incap;
	EREAL alpha_i = alpha;
	EREAL beta_i = beta;

	const EREAL *ap_i = ap;
	const EREAL *x_i = x;
	EREAL *y_i = y;
	EREAL head_rowsum, tail_rowsum;
	EREAL head_rowtmp, tail_rowtmp;
	EREAL matval;
	EREAL vecval;
	EREAL resval;
	EREAL head_tmp1, tail_tmp1;
	EREAL head_tmp2, tail_tmp2;
	FPU_FIX_DECL;

	incap = 1;




	if (incx < 0)
	  x_start = (-n + 1) * incx;
	else
	  x_start = 0;
	if (incy < 0)
	  y_start = (-n + 1) * incy;
	else
	  y_start = 0;

	if (n < 1) {
	  return;
	}
	if (alpha_i == 0.0 && beta_i == 1.0) {
	  return;
	}

	/* Check for error conditions. */
	if (order != blas_colmajor && order != blas_rowmajor) {
	  BLAS_error(routine_name, -1, order, NULL);
	}
	if (uplo != blas_upper && uplo != blas_lower) {
	  BLAS_error(routine_name, -2, uplo, NULL);
	}
	if (incx == 0) {
	  BLAS_error(routine_name, -7, incx, NULL);
	}
	if (incy == 0) {
	  BLAS_error(routine_name, -10, incy, NULL);
	}

	FPU_FIX_START;

	if (alpha_i == 0.0) {
	  {
	    y_index = y_start;
	    for (matrix_row = 0; matrix_row < n; matrix_row++) {
	      resval = y_i[y_index];

	      {
		/* Compute double_double = EREAL * EREAL. */
		EREAL a1, a2, b1, b2, con;

		con = beta_i * split;
		a1 = con - beta_i;
		a1 = con - a1;
		a2 = beta_i - a1;
		con = resval * split;
		b1 = con - resval;
		b1 = con - b1;
		b2 = resval - b1;

		head_tmp2 = beta_i * resval;
		tail_tmp2 =
		  (((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) + a2 * b2;
	      }

	      y_i[y_index] = head_tmp2;

	      y_index += incy;
	    }
	  }
	} else {
	  if (uplo == blas_lower)
	    order = (order == blas_rowmajor) ? blas_colmajor : blas_rowmajor;
	  if (order == blas_rowmajor) {
	    if (alpha_i == 1.0) {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    head_tmp1 = head_rowsum;
		    tail_tmp1 = tail_rowsum;
		    y_i[y_index] = head_tmp1;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    head_tmp1 = head_rowsum;
		    tail_tmp1 = tail_rowsum;
		    {
		      /* Compute double_double = EREAL * EREAL. */
		      EREAL a1, a2, b1, b2, con;

		      con = beta_i * split;
		      a1 = con - beta_i;
		      a1 = con - a1;
		      a2 = beta_i - a1;
		      con = resval * split;
		      b1 = con - resval;
		      b1 = con - b1;
		      b2 = resval - b1;

		      head_tmp2 = beta_i * resval;
		      tail_tmp2 =
			(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) +
			a2 * b2;
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
		      head_tmp2 = t1 + t2;
		      tail_tmp2 = t2 - (head_tmp2 - t1);
		    }
		    y_i[y_index] = head_tmp2;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      }
	    } else {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    {
		      /* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
		      EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_rowsum * split;
		      a11 = con - head_rowsum;
		      a11 = con - a11;
		      a21 = head_rowsum - a11;
		      con = alpha_i * split;
		      b1 = con - alpha_i;
		      b1 = con - b1;
		      b2 = alpha_i - b1;

		      c11 = head_rowsum * alpha_i;
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_rowsum * alpha_i;
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_tmp1 = t1 + t2;
		      tail_tmp1 = t2 - (head_tmp1 - t1);
		    }
		    y_i[y_index] = head_tmp1;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (n - step - 1) * incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    {
		      /* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
		      EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_rowsum * split;
		      a11 = con - head_rowsum;
		      a11 = con - a11;
		      a21 = head_rowsum - a11;
		      con = alpha_i * split;
		      b1 = con - alpha_i;
		      b1 = con - b1;
		      b2 = alpha_i - b1;

		      c11 = head_rowsum * alpha_i;
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_rowsum * alpha_i;
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_tmp1 = t1 + t2;
		      tail_tmp1 = t2 - (head_tmp1 - t1);
		    }
		    {
		      /* Compute double_double = EREAL * EREAL. */
		      EREAL a1, a2, b1, b2, con;

		      con = beta_i * split;
		      a1 = con - beta_i;
		      a1 = con - a1;
		      a2 = beta_i - a1;
		      con = resval * split;
		      b1 = con - resval;
		      b1 = con - b1;
		      b2 = resval - b1;

		      head_tmp2 = beta_i * resval;
		      tail_tmp2 =
			(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) +
			a2 * b2;
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
		      head_tmp2 = t1 + t2;
		      tail_tmp2 = t2 - (head_tmp2 - t1);
		    }
		    y_i[y_index] = head_tmp2;

		    y_index += incy;
		    ap_start += incap;
		  }
		}
	      }
	    }
	  } else {
	    if (alpha_i == 1.0) {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    head_tmp1 = head_rowsum;
		    tail_tmp1 = tail_rowsum;
		    y_i[y_index] = head_tmp1;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    head_tmp1 = head_rowsum;
		    tail_tmp1 = tail_rowsum;
		    {
		      /* Compute double_double = EREAL * EREAL. */
		      EREAL a1, a2, b1, b2, con;

		      con = beta_i * split;
		      a1 = con - beta_i;
		      a1 = con - a1;
		      a2 = beta_i - a1;
		      con = resval * split;
		      b1 = con - resval;
		      b1 = con - b1;
		      b2 = resval - b1;

		      head_tmp2 = beta_i * resval;
		      tail_tmp2 =
			(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) +
			a2 * b2;
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
		      head_tmp2 = t1 + t2;
		      tail_tmp2 = t2 - (head_tmp2 - t1);
		    }
		    y_i[y_index] = head_tmp2;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      }
	    } else {
	      if (beta_i == 0.0) {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    {
		      /* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
		      EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_rowsum * split;
		      a11 = con - head_rowsum;
		      a11 = con - a11;
		      a21 = head_rowsum - a11;
		      con = alpha_i * split;
		      b1 = con - alpha_i;
		      b1 = con - b1;
		      b2 = alpha_i - b1;

		      c11 = head_rowsum * alpha_i;
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_rowsum * alpha_i;
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_tmp1 = t1 + t2;
		      tail_tmp1 = t2 - (head_tmp1 - t1);
		    }
		    y_i[y_index] = head_tmp1;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      } else {
		{
		  y_index = y_start;
		  ap_start = 0;
		  for (matrix_row = 0; matrix_row < n; matrix_row++) {
		    x_index = x_start;
		    ap_index = ap_start;
		    head_rowsum = tail_rowsum = 0.0;
		    head_rowtmp = tail_rowtmp = 0.0;
		    for (step = 0; step < matrix_row; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += incap;
		      x_index += incx;
		    }
		    for (step = matrix_row; step < n; step++) {
		      matval = ap_i[ap_index];
		      vecval = x_i[x_index];
		      {
			/* Compute double_double = EREAL * EREAL. */
			EREAL a1, a2, b1, b2, con;

			con = matval * split;
			a1 = con - matval;
			a1 = con - a1;
			a2 = matval - a1;
			con = vecval * split;
			b1 = con - vecval;
			b1 = con - b1;
			b2 = vecval - b1;

			head_rowtmp = matval * vecval;
			tail_rowtmp =
			  (((a1 * b1 - head_rowtmp) + a1 * b2) + a2 * b1) +
			  a2 * b2;
		      }
		      {
			/* Compute EREAL-EREAL = EREAL-EREAL + EREAL-EREAL. */
			EREAL bv;
			EREAL s1, s2, t1, t2;

			/* Add two hi words. */
			s1 = head_rowsum + head_rowtmp;
			bv = s1 - head_rowsum;
			s2 = ((head_rowtmp - bv) + (head_rowsum - (s1 - bv)));

			/* Add two lo words. */
			t1 = tail_rowsum + tail_rowtmp;
			bv = t1 - tail_rowsum;
			t2 = ((tail_rowtmp - bv) + (tail_rowsum - (t1 - bv)));

			s2 += t1;

			/* Renormalize (s1, s2)  to  (t1, s2) */
			t1 = s1 + s2;
			s2 = s2 - (t1 - s1);

			t2 += s2;

			/* Renormalize (t1, t2)  */
			head_rowsum = t1 + t2;
			tail_rowsum = t2 - (head_rowsum - t1);
		      }
		      ap_index += (step + 1) * incap;
		      x_index += incx;
		    }
		    resval = y_i[y_index];
		    {
		      /* Compute EREAL-EREAL = EREAL-EREAL * EREAL. */
		      EREAL a11, a21, b1, b2, c11, c21, c2, con, t1, t2;

		      con = head_rowsum * split;
		      a11 = con - head_rowsum;
		      a11 = con - a11;
		      a21 = head_rowsum - a11;
		      con = alpha_i * split;
		      b1 = con - alpha_i;
		      b1 = con - b1;
		      b2 = alpha_i - b1;

		      c11 = head_rowsum * alpha_i;
		      c21 =
			(((a11 * b1 - c11) + a11 * b2) + a21 * b1) + a21 * b2;

		      c2 = tail_rowsum * alpha_i;
		      t1 = c11 + c2;
		      t2 = (c2 - (t1 - c11)) + c21;

		      head_tmp1 = t1 + t2;
		      tail_tmp1 = t2 - (head_tmp1 - t1);
		    }
		    {
		      /* Compute double_double = EREAL * EREAL. */
		      EREAL a1, a2, b1, b2, con;

		      con = beta_i * split;
		      a1 = con - beta_i;
		      a1 = con - a1;
		      a2 = beta_i - a1;
		      con = resval * split;
		      b1 = con - resval;
		      b1 = con - b1;
		      b2 = resval - b1;

		      head_tmp2 = beta_i * resval;
		      tail_tmp2 =
			(((a1 * b1 - head_tmp2) + a1 * b2) + a2 * b1) +
			a2 * b2;
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
		      head_tmp2 = t1 + t2;
		      tail_tmp2 = t2 - (head_tmp2 - t1);
		    }
		    y_i[y_index] = head_tmp2;

		    y_index += incy;
		    ap_start += (matrix_row + 1) * incap;
		  }
		}
	      }
	    }
	  }			/* if order == ... */
	}			/* alpha != 0 */

	FPU_FIX_STOP;
      }
      break;
    }

  }
}
