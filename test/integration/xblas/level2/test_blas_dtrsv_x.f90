program test_blas_dtrsv_x
    !! Triangular solve: x := alpha * op(T)^-1 * x.
    !! Use a diagonal-dominant T to keep the system well-conditioned.
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dtrsv_x, &
                              blas_upper, blas_no_trans, blas_non_unit_diag
    use ref_quad_xblas, only: ref_blas_dtrsv_x
    implicit none
    integer, parameter :: n = 64
    real(ep), allocatable :: t(:,:), x(:), x_ref(:), x_got(:)
    real(ep) :: alpha, err, tol
    integer :: i
    call report_init('blas_dtrsv_x', target_name)
    call gen_matrix_quad(n, n, t, seed = 11)
    do i = 1, n
        t(i, i) = t(i, i) + real(2 * n, ep)   ! diagonal dominance
    end do
    call gen_vector_quad(n, x, seed = 12)
    allocate(x_ref(n), x_got(n)); x_ref = x; x_got = x
    alpha = 1.7_ep
    call ref_blas_dtrsv_x(blas_upper, blas_no_trans, blas_non_unit_diag, &
                          n, alpha, t, n, x_ref, 1)
    call target_blas_dtrsv_x(blas_upper, blas_no_trans, blas_non_unit_diag, &
                             n, alpha, t, n, x_got, 1)
    err = max_rel_err_vec(x_got, x_ref)
    tol = 16.0_ep * (2.0_ep * real(n, ep) * real(n, ep)) * target_eps
    call report_case('n=64', err, tol)
    call report_finalize()
end program
