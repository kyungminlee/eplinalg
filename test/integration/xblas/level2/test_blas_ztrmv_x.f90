program test_blas_ztrmv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_matrix_complex, gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_ztrmv_x, &
                              blas_upper, blas_no_trans, blas_non_unit_diag
    use ref_quad_xblas, only: ref_blas_ztrmv_x
    implicit none
    integer, parameter :: n = 64
    complex(ep), allocatable :: t(:,:), x(:), x_ref(:), x_got(:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    call report_init('blas_ztrmv_x', target_name)
    call gen_matrix_complex(n, n, t, seed = 11)
    call gen_vector_complex(n, x, seed = 12)
    allocate(x_ref(n), x_got(n)); x_ref = x; x_got = x
    alpha = (1.7_ep, -0.4_ep)
    call ref_blas_ztrmv_x(blas_upper, blas_no_trans, blas_non_unit_diag, &
                          n, alpha, t, n, x_ref, 1)
    call target_blas_ztrmv_x(blas_upper, blas_no_trans, blas_non_unit_diag, &
                             n, alpha, t, n, x_got, 1)
    err = max_rel_err_vec_z(x_got, x_ref)
    tol = 16.0_ep * (4.0_ep * real(n, ep) * real(n, ep)) * target_eps
    call report_case('n=64', err, tol)
    call report_finalize()
end program
