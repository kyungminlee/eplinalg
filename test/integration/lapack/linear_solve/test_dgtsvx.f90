program test_dgtsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_vector_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgtsvx
    use ref_quad_lapack, only: dgtsvx
    implicit none
    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    real(ep), allocatable :: dl(:), d(:), du(:), B0(:,:)
    real(ep), allocatable :: dlf_r(:), df_r(:), duf_r(:), du2_r(:), X_r(:,:)
    real(ep), allocatable :: dlf_g(:), df_g(:), duf_g(:), du2_g(:), X_g(:,:)
    real(ep), allocatable :: ferr(:), berr(:), work(:)
    integer,  allocatable :: ipiv_r(:), ipiv_g(:), iwork(:)
    real(ep) :: rcond_r, rcond_g, err, tol
    character(len=48) :: label
    call report_init('dgtsvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_quad(n,   d,  seed = 90001 + 19 * i)
        call gen_vector_quad(n-1, dl, seed = 90011 + 19 * i)
        call gen_vector_quad(n-1, du, seed = 90021 + 19 * i)
        do j = 1, n; d(j) = d(j) + 4.0_ep; end do
        call gen_matrix_quad(n, nrhs, B0, seed = 90031 + 19 * i)
        allocate(dlf_r(max(1,n-1)), df_r(n), duf_r(max(1,n-1)), du2_r(max(1,n-2)), X_r(n, nrhs))
        allocate(dlf_g(max(1,n-1)), df_g(n), duf_g(max(1,n-1)), du2_g(max(1,n-2)), X_g(n, nrhs))
        allocate(ipiv_r(n), ipiv_g(n), iwork(n), work(3*n), ferr(nrhs), berr(nrhs))
        call dgtsvx('N', 'N', n, nrhs, dl, d, du, dlf_r, df_r, duf_r, du2_r, ipiv_r, &
                    B0, n, X_r, n, rcond_r, ferr, berr, work, iwork, info)
        call target_dgtsvx('N', 'N', n, nrhs, dl, d, du, dlf_g, df_g, duf_g, du2_g, ipiv_g, &
                           B0, n, X_g, n, rcond_g, ferr, berr, info)
        err = max(max_rel_err_mat(X_g, X_r), rel_err_scalar(rcond_g, rcond_r))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dlf_r, df_r, duf_r, du2_r, X_r, dlf_g, df_g, duf_g, du2_g, X_g, ipiv_r, ipiv_g, iwork, work, ferr, berr)
    end do
    call report_finalize()
end program test_dgtsvx
