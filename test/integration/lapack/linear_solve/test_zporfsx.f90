! zporfsx: extra-precise refinement of HPD Cholesky solution.
program test_zporfsx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zpotrf, target_zpotrs, target_zporfsx
    use ref_quad_lapack, only: zpotrf, zpotrs, zporfsx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, info
    complex(ep), allocatable :: A0(:,:), B0(:,:)
    complex(ep), allocatable :: A_ref(:,:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: A_got(:,:), B_got(:,:), X_got(:,:), work(:)
    real(ep),    allocatable :: S(:), berr_ref(:), berr_got(:),               &
                                ebn_ref(:,:), ebn_got(:,:), ebc_ref(:,:), ebc_got(:,:), &
                                rwork(:), params(:)
    real(ep)    :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zporfsx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 411101 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 411111 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,nrhs), X_ref(n,nrhs))
        allocate(A_got(n,n), B_got(n,nrhs), X_got(n,nrhs))
        allocate(S(n), berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(work(2*n), rwork(2*n), params(1))
        A_ref = A0; A_got = A0
        call zpotrf('U', n, A_ref, n, info)
        call target_zpotrf('U', n, A_got, n, info)
        X_ref = B0; X_got = B0
        call zpotrs('U', n, nrhs, A_ref, n, X_ref, n, info)
        call target_zpotrs('U', n, nrhs, A_got, n, X_got, n, info)
        S = 0.0_ep
        call zporfsx('U', 'N', n, nrhs, A0, n, A_ref, n, S,                &
                     B0, n, X_ref, n, rcond_ref, berr_ref,                 &
                     n_err_bnds, ebn_ref, ebc_ref, 0, params, work, rwork, info)
        call target_zporfsx('U', 'N', n, nrhs, A0, n, A_got, n, S,         &
                            B0, n, X_got, n, rcond_got, berr_got,          &
                            n_err_bnds, ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat_z(X_got, X_ref),                           &
                  rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, X_ref, A_got, B_got, X_got, S,    &
                   berr_ref, berr_got, ebn_ref, ebn_got, ebc_ref, ebc_got, &
                   work, rwork, params)
    end do
    call report_finalize()
end program test_zporfsx
