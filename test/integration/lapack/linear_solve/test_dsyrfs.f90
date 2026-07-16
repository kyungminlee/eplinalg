program test_dsyrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dsytrf, target_dsytrs, target_dsyrfs
    use ref_quad_lapack, only: dsytrf, dsytrs, dsyrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, lwork
    integer, allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), B0(:,:), X_ref(:,:), X_got(:,:)
    real(ep), allocatable :: work_f(:), work(:), ferr(:), berr(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsyrfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 146001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 146011 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(3*n), iwork(n), ferr(nrhs), berr(nrhs))
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
        call dsyrfs('U', n, nrhs, A0, n, A_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, ferr, berr, work, iwork, info)
        call target_dsyrfs('U', n, nrhs, A0, n, A_got, n, ipiv_got, B0, n, &
                           X_got, n, ferr, berr, info)
        err = max_rel_err_mat(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, X_ref, X_got, work, iwork, ferr, berr)
    end do
    call report_finalize()
end program test_dsyrfs
