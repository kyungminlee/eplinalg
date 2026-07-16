! zgbrfsx: extra-precise refinement of complex banded LU solution.
program test_zgbrfsx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zgbtrf, target_zgbtrs, target_zgbrfsx
    use ref_quad_lapack, only: zgbtrf, zgbtrs, zgbrfsx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: kl           = 2
    integer, parameter :: ku           = 3
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, ldab, ldafb, info, j, k
    complex(ep), allocatable :: Adense(:,:), AB0(:,:), AFB_ref(:,:), AFB_got(:,:)
    complex(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep),    allocatable :: R(:), C(:), berr_ref(:), berr_got(:),    &
                                ebn_ref(:,:), ebn_got(:,:), ebc_ref(:,:), ebc_got(:,:), &
                                rwork(:), params(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep)    :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zgbrfsx', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kl + ku + 1; ldafb = 2*kl + ku + 1
        call gen_matrix_complex(n, n, Adense, seed = 410701 + 47 * i)
        do j = 1, n; Adense(j, j) = Adense(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        allocate(AB0(ldab, n)); AB0 = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        allocate(AFB_ref(ldafb, n), AFB_got(ldafb, n)); AFB_ref = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AFB_ref(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        AFB_got = AFB_ref
        call gen_matrix_complex(n, nrhs, B0, seed = 410711 + 47 * i)
        allocate(ipiv_ref(n), ipiv_got(n), X_ref(n, nrhs), X_got(n, nrhs))
        allocate(R(n), C(n), berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(work(2*n), rwork(2*n), params(1))
        call zgbtrf(n, n, kl, ku, AFB_ref, ldafb, ipiv_ref, info)
        call target_zgbtrf(n, n, kl, ku, AFB_got, ldafb, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zgbtrs('N', n, kl, ku, nrhs, AFB_ref, ldafb, ipiv_ref, X_ref, n, info)
        call target_zgbtrs('N', n, kl, ku, nrhs, AFB_got, ldafb, ipiv_got, X_got, n, info)
        R = 0.0_ep; C = 0.0_ep
        call zgbrfsx('N', 'N', n, kl, ku, nrhs, AB0, ldab, AFB_ref, ldafb,        &
                     ipiv_ref, R, C, B0, n, X_ref, n, rcond_ref, berr_ref,        &
                     n_err_bnds, ebn_ref, ebc_ref, 0, params, work, rwork, info)
        call target_zgbrfsx('N', 'N', n, kl, ku, nrhs, AB0, ldab, AFB_got, ldafb, &
                            ipiv_got, R, C, B0, n, X_got, n, rcond_got, berr_got, &
                            n_err_bnds, ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat_z(X_got, X_ref),                                  &
                  rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AFB_ref, AFB_got, B0, ipiv_ref, ipiv_got, &
                   X_ref, X_got, R, C, berr_ref, berr_got,                &
                   ebn_ref, ebn_got, ebc_ref, ebc_got, work, rwork, params)
    end do
    call report_finalize()
end program test_zgbrfsx
