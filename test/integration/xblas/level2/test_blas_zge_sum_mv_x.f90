program test_blas_zge_sum_mv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_matrix_complex, gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zge_sum_mv_x
    use ref_quad_xblas, only: ref_blas_zge_sum_mv_x
    implicit none
    integer, parameter :: m = 64, n = 48
    complex(ep), allocatable :: a(:,:), b(:,:), x(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    call report_init('blas_zge_sum_mv_x', target_name)
    call gen_matrix_complex(m, n, a, seed = 11)
    call gen_matrix_complex(m, n, b, seed = 12)
    call gen_vector_complex(n, x, seed = 13)
    call gen_vector_complex(m, y_ref, seed = 14)
    allocate(y_got(m)); y_got = y_ref
    alpha = (1.7_ep, -0.4_ep); beta = (-0.3_ep, 0.6_ep)
    call ref_blas_zge_sum_mv_x(m, n, alpha, a, m, x, 1, beta, b, m, y_ref, 1)
    call target_blas_zge_sum_mv_x(m, n, alpha, a, m, x, 1, beta, b, m, y_got, 1)
    err = max_rel_err_vec_z(y_got, y_ref)
    tol = 16.0_ep * (8.0_ep * real(m, ep) * real(n, ep)) * target_eps
    call report_case('m=64 n=48', err, tol)
    call report_finalize()
end program
