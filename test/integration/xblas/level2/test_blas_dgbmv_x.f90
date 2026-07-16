program test_blas_dgbmv_x
    !! Banded GEMV: A is m×n with kl sub-diagonals and ku super-diagonals.
    !! Storage: column j is packed into A_band(:, j), with diagonal at
    !! row ku+1 (BLAS convention).
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dgbmv_x, &
                              blas_no_trans
    use ref_quad_xblas, only: ref_blas_dgbmv_x
    implicit none

    integer, parameter :: m = 64, n = 64, kl = 3, ku = 2
    integer, parameter :: lda = kl + ku + 1
    real(ep), allocatable :: a(:,:), x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('blas_dgbmv_x', target_name)

    call gen_matrix_quad(lda, n, a, seed = 1011)
    call gen_vector_quad(n, x,    seed = 1012)
    call gen_vector_quad(m, y_ref, seed = 1013)
    allocate(y_got(m))
    y_got = y_ref
    alpha = 1.7_ep
    beta  = -0.3_ep

    call ref_blas_dgbmv_x(blas_no_trans, m, n, kl, ku, alpha, a, lda, &
                          x, 1, beta, y_ref, 1)
    call target_blas_dgbmv_x(blas_no_trans, m, n, kl, ku, alpha, a, lda, &
                             x, 1, beta, y_got, 1)

    err = max_rel_err_vec(y_got, y_ref)
    tol = 16.0_ep * (2.0_ep * real(m, ep) * real(kl + ku + 1, ep)) * target_eps
    write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm=', m, ' n=', n, ' kl=', kl, ' ku=', ku
    call report_case(trim(label), err, tol)

    call report_finalize()
end program test_blas_dgbmv_x
