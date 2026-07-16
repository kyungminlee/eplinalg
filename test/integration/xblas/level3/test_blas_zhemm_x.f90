program test_blas_zhemm_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_mat_z
    use test_data,      only: gen_matrix_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zhemm_x, &
                              blas_left_side, blas_upper
    use ref_quad_xblas, only: ref_blas_zhemm_x
    implicit none
    integer, parameter :: m = 32, n = 24
    complex(ep), allocatable :: a(:,:), b(:,:), c_ref(:,:), c_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    integer :: i
    call report_init('blas_zhemm_x', target_name)
    call gen_matrix_complex(m, m, a, seed = 11)
    do i = 1, m
        a(i, i) = cmplx(real(a(i, i), ep), 0.0_ep, ep)
    end do
    call gen_matrix_complex(m, n, b, seed = 12)
    call gen_matrix_complex(m, n, c_ref, seed = 13)
    allocate(c_got(m, n)); c_got = c_ref
    alpha = (1.7_ep, -0.4_ep); beta = (-0.3_ep, 0.6_ep)
    call ref_blas_zhemm_x(blas_left_side, blas_upper, m, n, alpha, a, m, &
                          b, m, beta, c_ref, m)
    call target_blas_zhemm_x(blas_left_side, blas_upper, m, n, alpha, a, m, &
                             b, m, beta, c_got, m)
    err = max_rel_err_mat_z(c_got, c_ref)
    tol = 16.0_ep * (4.0_ep * real(m, ep) * real(m, ep) * real(n, ep)) * target_eps
    call report_case('m=32 n=24', err, tol)
    call report_finalize()
end program
