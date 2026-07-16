program test_zposvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zposvx
    use ref_quad_lapack, only: zposvx
    implicit none

    integer, parameter :: ns(*)   = [16, 32, 48]
    integer, parameter :: nrhs    = 2
    integer :: i, n, info
    complex(ep), allocatable :: A0(:,:), B0(:,:)
    complex(ep), allocatable :: A_ref(:,:), AF_ref(:,:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: A_got(:,:), AF_got(:,:), B_got(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: S_ref(:), S_got(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), rwork(:)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('zposvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 400061 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 400071 + 47 * i)
        allocate(A_ref(n,n), AF_ref(n,n), B_ref(n,nrhs), X_ref(n,nrhs), S_ref(n))
        allocate(A_got(n,n), AF_got(n,n), B_got(n,nrhs), X_got(n,nrhs), S_got(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(work(2*n), rwork(n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call zposvx('E', 'U', n, nrhs, A_ref, n, AF_ref, n, equed_ref, S_ref, &
                    B_ref, n, X_ref, n, rcond_ref, ferr_ref, berr_ref, work, rwork, info)
        call target_zposvx('E', 'U', n, nrhs, A_got, n, AF_got, n, equed_got, S_got, &
                           B_got, n, X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat_z(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, AF_ref, B_ref, X_ref, S_ref, &
                   A_got, AF_got, B_got, X_got, S_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, work, rwork)
    end do
    call report_finalize()
end program test_zposvx
