program test_zhpsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hermitian_matrix_quad, gen_matrix_complex, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhpsvx
    use ref_quad_lapack, only: zhpsvx
    implicit none
    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: nrhs  = 2
    integer :: i, n, np, info, j
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), AP(:), AFP_ref(:), AFP_got(:), B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep),    allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label
    call report_init('zhpsvx', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A0, seed = 76001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 76011 + 47 * i)
        allocate(AP(np), AFP_ref(np), AFP_got(np), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        call pack_herm_packed_quad('U', n, A0, AP)
        call zhpsvx('N', 'U', n, nrhs, AP, AFP_ref, ipiv_ref, B0, n, X_ref, n, &
                    rcond_ref, ferr, berr, work, rwork, info)
        call target_zhpsvx('N', 'U', n, nrhs, AP, AFP_got, ipiv_got, B0, n, X_got, n, &
                           rcond_got, ferr, berr, info)
        err = max(max_rel_err_mat_z(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP, AFP_ref, AFP_got, ipiv_ref, ipiv_got, X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zhpsvx
