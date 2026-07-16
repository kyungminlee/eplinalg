program test_blas_dgemv2_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dgemv2_x, &
                              blas_no_trans
    use ref_quad_xblas, only: ref_blas_dgemv2_x
    implicit none
    integer, parameter :: m = 64, n = 48
    real(ep), allocatable :: a(:,:), head_x(:), tail_x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    call report_init('blas_dgemv2_x', target_name)
    call gen_matrix_quad(m, n, a, seed = 11)
    call gen_vector_quad(n, head_x, seed = 12)
    call gen_vector_quad(n, tail_x, seed = 13)
    call gen_vector_quad(m, y_ref, seed = 14)
    allocate(y_got(m)); y_got = y_ref
    alpha = 1.7_ep; beta = -0.3_ep
    call ref_blas_dgemv2_x(blas_no_trans, m, n, alpha, a, m, head_x, tail_x, 1, &
                           beta, y_ref, 1)
    call target_blas_dgemv2_x(blas_no_trans, m, n, alpha, a, m, head_x, tail_x, 1, &
                              beta, y_got, 1)
    err = max_rel_err_vec(y_got, y_ref)
    tol = 16.0_ep * (4.0_ep * real(m, ep) * real(n, ep)) * target_eps
    call report_case('m=64 n=48', err, tol)
    call report_finalize()
end program
