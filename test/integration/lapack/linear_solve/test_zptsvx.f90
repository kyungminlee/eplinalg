program test_zptsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_vector_quad, gen_vector_complex, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zptsvx
    use ref_quad_lapack, only: zptsvx
    implicit none
    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    real(ep),    allocatable :: d(:), df_r(:), df_g(:), rwork(:), ferr(:), berr(:)
    complex(ep), allocatable :: e(:), B0(:,:), ef_r(:), ef_g(:), X_r(:,:), X_g(:,:), work(:)
    real(ep) :: rcond_r, rcond_g, err, tol
    character(len=48) :: label
    call report_init('zptsvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n, d, seed = 93001 + 19 * i)
        call gen_vector_complex(n-1, e, seed = 93011 + 19 * i)
        do j = 1, n; d(j) = abs(d(j)) + 4.0_ep; end do
        call gen_matrix_complex(n, nrhs, B0, seed = 93021 + 19 * i)
        allocate(df_r(n), ef_r(max(1,n-1)), X_r(n, nrhs))
        allocate(df_g(n), ef_g(max(1,n-1)), X_g(n, nrhs))
        allocate(work(n), rwork(n), ferr(nrhs), berr(nrhs))
        call zptsvx('N', n, nrhs, d, e, df_r, ef_r, B0, n, X_r, n, rcond_r, ferr, berr, work, rwork, info)
        call target_zptsvx('N', n, nrhs, d, e, df_g, ef_g, B0, n, X_g, n, rcond_g, ferr, berr, info)
        err = max(max_rel_err_mat_z(X_g, X_r), rel_err_scalar(rcond_g, rcond_r))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(df_r, ef_r, X_r, df_g, ef_g, X_g, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zptsvx
