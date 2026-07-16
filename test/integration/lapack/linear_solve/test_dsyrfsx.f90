! dsyrfsx: extra-precise refinement of symmetric Bunch–Kaufman solution.
program test_dsyrfsx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dsytrf, target_dsytrs, target_dsyrfsx
    use ref_quad_lapack, only: dsytrf, dsytrs, dsyrfsx
    implicit none

    integer, parameter :: ns(*)        = [16, 32, 48]
    integer, parameter :: nrhs         = 2
    integer, parameter :: n_err_bnds   = 3
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), B0(:,:), X_ref(:,:), X_got(:,:)
    real(ep), allocatable :: work_f(:), S(:), berr_ref(:), berr_got(:),       &
                              ebn_ref(:,:), ebn_got(:,:), ebc_ref(:,:), ebc_got(:,:), &
                              work(:), params(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep)  :: rcond_ref, rcond_got, wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsyrfsx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 411401 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 411411 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs))
        allocate(S(n), berr_ref(nrhs), berr_got(nrhs))
        allocate(ebn_ref(nrhs, n_err_bnds), ebn_got(nrhs, n_err_bnds))
        allocate(ebc_ref(nrhs, n_err_bnds), ebc_got(nrhs, n_err_bnds))
        allocate(work(4*n), iwork(n), params(1))
        A_ref = A0; A_got = A0
        call dsytrf('U', n, A_ref, n, ipiv_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work_f(lwork))
        call dsytrf('U', n, A_ref, n, ipiv_ref, work_f, lwork, info)
        deallocate(work_f)
        call target_dsytrf('U', n, A_got, n, ipiv_got, info)
        X_ref = B0; X_got = B0
        call dsytrs('U', n, nrhs, A_ref, n, ipiv_ref, X_ref, n, info)
        call target_dsytrs('U', n, nrhs, A_got, n, ipiv_got, X_got, n, info)
        S = 0.0_ep
        call dsyrfsx('U', 'N', n, nrhs, A0, n, A_ref, n, ipiv_ref, S,       &
                     B0, n, X_ref, n, rcond_ref, berr_ref,                  &
                     n_err_bnds, ebn_ref, ebc_ref, 0, params, work, iwork, info)
        call target_dsyrfsx('U', 'N', n, nrhs, A0, n, A_got, n, ipiv_got, S, &
                            B0, n, X_got, n, rcond_got, berr_got,            &
                            n_err_bnds, ebn_got, ebc_got, 0, params, info)
        err = max(max_rel_err_mat(X_got, X_ref),                            &
                  rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, X_ref, X_got, S,       &
                   berr_ref, berr_got, ebn_ref, ebn_got, ebc_ref, ebc_got,  &
                   work, iwork, params)
    end do
    call report_finalize()
end program test_dsyrfsx
