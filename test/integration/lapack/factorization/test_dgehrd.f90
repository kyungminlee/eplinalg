program test_dgehrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgehrd
    use ref_quad_lapack, only: dgehrd
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    real(ep) :: wopt(1), err_a, err_t, tol
    character(len=48) :: label

    call report_init('dgehrd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 36001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), tau_ref(n-1), tau_got(n-1))
        A_ref = A0; A_got = A0
        call dgehrd(n, 1, n, A_ref, n, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A_ref, n, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgehrd(n, 1, n, A_got, n, tau_got, info)
        err_a = max_rel_err_mat(A_got, A_ref)
        err_t = max_rel_err_vec(tau_got, tau_ref)
        tol   = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a)') 'n=', n, ',out=A'
        call report_case(trim(label), err_a, tol)
        write(label, '(a,i0,a)') 'n=', n, ',out=tau'
        call report_case(trim(label), err_t, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dgehrd
