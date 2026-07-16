! zgbsvxx: extra-precise expert solver, complex banded.
program test_zgbsvxx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgbsvxx
    use ref_quad_lapack, only: zgbsvxx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: kl           = 2
    integer, parameter :: ku           = 3
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, ldab, ldafb, info, j, k
    complex(ep), allocatable :: A(:,:), B0(:,:)
    complex(ep), allocatable :: AB_ref(:,:), AFB_ref(:,:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: AB_got(:,:), AFB_got(:,:), B_got(:,:), X_got(:,:), work(:)
    real(ep),    allocatable :: R_ref(:), C_ref(:), R_got(:), C_got(:)
    real(ep),    allocatable :: berr_ref(:), berr_got(:), ebn_ref(:,:), ebn_got(:,:), &
                                ebc_ref(:,:), ebc_got(:,:), rwork(:), params(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    real(ep)    :: rcond_ref, rcond_got, rpvgrw_ref, rpvgrw_got, err, tol
    character   :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('zgbsvxx', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kl + ku + 1; ldafb = 2*kl + ku + 1
        call gen_matrix_complex(n, n, A, seed = 410501 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 410511 + 47 * i)
        allocate(AB_ref(ldab, n), AFB_ref(ldafb, n), B_ref(n, nrhs), X_ref(n, nrhs), R_ref(n), C_ref(n))
        allocate(AB_got(ldab, n), AFB_got(ldafb, n), B_got(n, nrhs), X_got(n, nrhs), R_got(n), C_got(n))
        allocate(berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(ipiv_ref(n), ipiv_got(n), work(2*n), rwork(2*n), params(1))
        AB_ref = (0.0_ep, 0.0_ep)
        do j = 1, n
            do k = max(1, j - ku), min(n, j + kl)
                AB_ref(ku + 1 + k - j, j) = A(k, j)
            end do
        end do
        AB_got = AB_ref
        B_ref = B0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call zgbsvxx('E', 'N', n, kl, ku, nrhs, AB_ref, ldab, AFB_ref, ldafb, &
                     ipiv_ref, equed_ref, R_ref, C_ref, B_ref, n, X_ref, n,    &
                     rcond_ref, rpvgrw_ref, berr_ref, n_err_bnds,              &
                     ebn_ref, ebc_ref, 0, params, work, rwork, info)
        call target_zgbsvxx('E', 'N', n, kl, ku, nrhs, AB_got, ldab, AFB_got, ldafb, &
                            ipiv_got, equed_got, R_got, C_got, B_got, n, X_got, n,   &
                            rcond_got, rpvgrw_got, berr_got, n_err_bnds,             &
                            ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat_z(X_got, X_ref),         &
                  rel_err_scalar(rcond_got, rcond_ref),  &
                  rel_err_scalar(rpvgrw_got, rpvgrw_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AB_ref, AFB_ref, B_ref, X_ref, R_ref, C_ref, &
                   AB_got, AFB_got, B_got, X_got, R_got, C_got,        &
                   berr_ref, berr_got, ebn_ref, ebn_got, ebc_ref, ebc_got, &
                   ipiv_ref, ipiv_got, work, rwork, params)
    end do
    call report_finalize()
end program test_zgbsvxx
