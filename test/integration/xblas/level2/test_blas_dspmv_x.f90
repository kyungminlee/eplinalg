program test_blas_dspmv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dspmv_x, &
                              blas_upper
    use ref_quad_xblas, only: ref_blas_dspmv_x
    implicit none
    integer, parameter :: n = 64
    integer, parameter :: lap = n * (n + 1) / 2
    real(ep), allocatable :: ap(:), x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    call report_init('blas_dspmv_x', target_name)
    call gen_vector_quad(lap, ap, seed = 11)
    call gen_vector_quad(n, x, seed = 12)
    call gen_vector_quad(n, y_ref, seed = 13)
    allocate(y_got(n)); y_got = y_ref
    alpha = 1.7_ep; beta = -0.3_ep
    call ref_blas_dspmv_x(blas_upper, n, alpha, ap, x, 1, beta, y_ref, 1)
    call target_blas_dspmv_x(blas_upper, n, alpha, ap, x, 1, beta, y_got, 1)
    err = max_rel_err_vec(y_got, y_ref)
    tol = 16.0_ep * (2.0_ep * real(n, ep) * real(n, ep)) * target_eps
    call report_case('n=64', err, tol)
    call report_finalize()
end program
