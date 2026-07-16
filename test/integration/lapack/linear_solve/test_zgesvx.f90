program test_zgesvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgesvx
    use ref_quad_lapack, only: zgesvx
    implicit none

    integer, parameter :: ns(*)   = [16, 32, 48]
    integer, parameter :: nrhs    = 2
    integer :: i, n, j, info
    complex(ep), allocatable :: A0(:,:), B0(:,:)
    complex(ep), allocatable :: A_ref(:,:), AF_ref(:,:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: A_got(:,:), AF_got(:,:), B_got(:,:), X_got(:,:)
    complex(ep), allocatable :: work(:)
    real(ep), allocatable :: R_ref(:), C_ref(:), R_got(:), C_got(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), rwork(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('zgesvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 400021 + 47 * i)
        do j = 1, n; A0(j, j) = A0(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 400031 + 47 * i)
        allocate(A_ref(n,n), AF_ref(n,n), B_ref(n,nrhs), X_ref(n,nrhs), R_ref(n), C_ref(n))
        allocate(A_got(n,n), AF_got(n,n), B_got(n,nrhs), X_got(n,nrhs), R_got(n), C_got(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(ipiv_ref(n), ipiv_got(n), work(2*n), rwork(2*n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call zgesvx('E', 'N', n, nrhs, A_ref, n, AF_ref, n, ipiv_ref, equed_ref, &
                    R_ref, C_ref, B_ref, n, X_ref, n, rcond_ref, ferr_ref, berr_ref, &
                    work, rwork, info)
        call target_zgesvx('E', 'N', n, nrhs, A_got, n, AF_got, n, ipiv_got, equed_got, &
                           R_got, C_got, B_got, n, X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat_z(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, AF_ref, B_ref, X_ref, R_ref, C_ref, &
                   A_got, AF_got, B_got, X_got, R_got, C_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, ipiv_ref, ipiv_got, work, rwork)
    end do
    call report_finalize()
end program test_zgesvx
