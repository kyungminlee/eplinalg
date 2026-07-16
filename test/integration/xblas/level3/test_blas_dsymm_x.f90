program test_blas_dsymm_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_mat
    use test_data,      only: gen_matrix_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dsymm_x, &
                              blas_left_side, blas_upper
    use ref_quad_xblas, only: ref_blas_dsymm_x
    implicit none
    integer, parameter :: m = 32, n = 24
    real(ep), allocatable :: a(:,:), b(:,:), c_ref(:,:), c_got(:,:)
    real(ep) :: alpha, beta, err, tol
    call report_init('blas_dsymm_x', target_name)
    call gen_matrix_quad(m, m, a, seed = 11)   ! left side: A is m×m symmetric
    call gen_matrix_quad(m, n, b, seed = 12)
    call gen_matrix_quad(m, n, c_ref, seed = 13)
    allocate(c_got(m, n)); c_got = c_ref
    alpha = 1.7_ep; beta = -0.3_ep
    call ref_blas_dsymm_x(blas_left_side, blas_upper, m, n, alpha, a, m, &
                          b, m, beta, c_ref, m)
    call target_blas_dsymm_x(blas_left_side, blas_upper, m, n, alpha, a, m, &
                             b, m, beta, c_got, m)
    err = max_rel_err_mat(c_got, c_ref)
    tol = 16.0_ep * (2.0_ep * real(m, ep) * real(m, ep) * real(n, ep)) * target_eps
    call report_case('m=32 n=24', err, tol)
    call report_finalize()
end program
