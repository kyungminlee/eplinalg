program test_blas_dtbsv_x
    !! Banded triangular solve. Use unit-diagonal so we don't have to
    !! patch the band storage for diagonal dominance.
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dtbsv_x, &
                              blas_upper, blas_no_trans, blas_unit_diag
    use ref_quad_xblas, only: ref_blas_dtbsv_x
    implicit none
    integer, parameter :: n = 64, k = 3, ldt = k + 1
    real(ep), allocatable :: t(:,:), x(:), x_ref(:), x_got(:)
    real(ep) :: alpha, err, tol
    integer :: j
    call report_init('blas_dtbsv_x', target_name)
    call gen_matrix_quad(ldt, n, t, seed = 11)
    ! Scale off-diagonal bands to be small so the unit-diagonal solve is stable.
    do j = 1, n
        t(1:k, j) = t(1:k, j) * 0.1_ep
    end do
    call gen_vector_quad(n, x, seed = 12)
    allocate(x_ref(n), x_got(n)); x_ref = x; x_got = x
    alpha = 1.7_ep
    call ref_blas_dtbsv_x(blas_upper, blas_no_trans, blas_unit_diag, &
                          n, k, alpha, t, ldt, x_ref, 1)
    call target_blas_dtbsv_x(blas_upper, blas_no_trans, blas_unit_diag, &
                             n, k, alpha, t, ldt, x_got, 1)
    err = max_rel_err_vec(x_got, x_ref)
    tol = 16.0_ep * (2.0_ep * real(n, ep) * real(k + 1, ep)) * target_eps
    call report_case('n=64 k=3', err, tol)
    call report_finalize()
end program
