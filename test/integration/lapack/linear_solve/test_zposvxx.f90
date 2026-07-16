! zposvxx: extra-precise expert solver for HPD systems.
program test_zposvxx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zposvxx
    use ref_quad_lapack, only: zposvxx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, info
    complex(ep), allocatable :: A0(:,:), B0(:,:)
    complex(ep), allocatable :: A_ref(:,:), AF_ref(:,:), B_ref(:,:), X_ref(:,:)
    complex(ep), allocatable :: A_got(:,:), AF_got(:,:), B_got(:,:), X_got(:,:), work(:)
    real(ep),    allocatable :: S_ref(:), S_got(:), berr_ref(:), berr_got(:), &
                                ebn_ref(:,:), ebn_got(:,:), ebc_ref(:,:), ebc_got(:,:), &
                                rwork(:), params(:)
    real(ep)    :: rcond_ref, rcond_got, rpvgrw_ref, rpvgrw_got, err, tol
    character   :: equed_ref, equed_got
    character(len=48) :: label

    call report_init('zposvxx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 410901 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 410911 + 47 * i)
        allocate(A_ref(n,n), AF_ref(n,n), B_ref(n,nrhs), X_ref(n,nrhs), S_ref(n))
        allocate(A_got(n,n), AF_got(n,n), B_got(n,nrhs), X_got(n,nrhs), S_got(n))
        allocate(berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(work(2*n), rwork(2*n), params(1))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        equed_ref = 'N'; equed_got = 'N'
        call zposvxx('E', 'U', n, nrhs, A_ref, n, AF_ref, n, equed_ref, S_ref, &
                     B_ref, n, X_ref, n, rcond_ref, rpvgrw_ref, berr_ref,      &
                     n_err_bnds, ebn_ref, ebc_ref, 0, params, work, rwork, info)
        call target_zposvxx('E', 'U', n, nrhs, A_got, n, AF_got, n, equed_got, S_got, &
                            B_got, n, X_got, n, rcond_got, rpvgrw_got, berr_got,      &
                            n_err_bnds, ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat_z(X_got, X_ref),         &
                  rel_err_scalar(rcond_got, rcond_ref),  &
                  rel_err_scalar(rpvgrw_got, rpvgrw_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, AF_ref, B_ref, X_ref, S_ref,         &
                   A_got, AF_got, B_got, X_got, S_got,                 &
                   berr_ref, berr_got, ebn_ref, ebn_got, ebc_ref, ebc_got, &
                   work, rwork, params)
    end do
    call report_finalize()
end program test_zposvxx
