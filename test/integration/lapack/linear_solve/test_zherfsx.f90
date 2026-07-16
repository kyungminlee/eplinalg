! zherfsx: extra-precise refinement of Hermitian Bunch–Kaufman solution.
program test_zherfsx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zhetrf, target_zhetrs, target_zherfsx
    use ref_quad_lapack, only: zhetrf, zhetrs, zherfsx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work_f(:), work(:)
    real(ep),    allocatable :: S(:), berr_ref(:), berr_got(:),               &
                                ebn_ref(:,:), ebn_got(:,:), ebc_ref(:,:), ebc_got(:,:), &
                                rwork(:), params(:)
    integer,     allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep)    :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zherfsx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 411701 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 411711 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs))
        allocate(S(n), berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(work(2*n), rwork(2*n), params(1))
        A_ref = A0; A_got = A0
        call zhetrf('U', n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work_f(lwork))
        call zhetrf('U', n, A_ref, n, ipiv_ref, work_f, lwork, info)
        deallocate(work_f)
        call target_zhetrf('U', n, A_got, n, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zhetrs('U', n, nrhs, A_ref, n, ipiv_ref, X_ref, n, info)
        call target_zhetrs('U', n, nrhs, A_got, n, ipiv_got, X_got, n, info)
        S = 0.0_ep
        call zherfsx('U', 'N', n, nrhs, A0, n, A_ref, n, ipiv_ref, S,       &
                     B0, n, X_ref, n, rcond_ref, berr_ref,                  &
                     n_err_bnds, ebn_ref, ebc_ref, 0, params, work, rwork, info)
        call target_zherfsx('U', 'N', n, nrhs, A0, n, A_got, n, ipiv_got, S, &
                            B0, n, X_got, n, rcond_got, berr_got,            &
                            n_err_bnds, ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat_z(X_got, X_ref),                            &
                  rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, X_ref, X_got, S,       &
                   berr_ref, berr_got, ebn_ref, ebn_got, ebc_ref, ebc_got,  &
                   work, rwork, params)
    end do
    call report_finalize()
end program test_zherfsx
