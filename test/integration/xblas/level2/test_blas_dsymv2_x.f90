program test_blas_dsymv2_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dsymv2_x, &
                              blas_upper
    use ref_quad_xblas, only: ref_blas_dsymv2_x
    implicit none
    integer, parameter :: n = 64
    real(ep), allocatable :: a(:,:), head_x(:), tail_x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    call report_init('blas_dsymv2_x', target_name)
    call gen_matrix_quad(n, n, a, seed = 11)
    call gen_vector_quad(n, head_x, seed = 12)
    call gen_vector_quad(n, tail_x, seed = 13)
    call gen_vector_quad(n, y_ref, seed = 14)
    allocate(y_got(n)); y_got = y_ref
    alpha = 1.7_ep; beta = -0.3_ep
    call ref_blas_dsymv2_x(blas_upper, n, alpha, a, n, head_x, tail_x, 1, &
                           beta, y_ref, 1)
    call target_blas_dsymv2_x(blas_upper, n, alpha, a, n, head_x, tail_x, 1, &
                              beta, y_got, 1)
    err = max_rel_err_vec(y_got, y_ref)
    tol = 16.0_ep * (4.0_ep * real(n, ep) * real(n, ep)) * target_eps
    call report_case('n=64', err, tol)
    call report_finalize()
end program
