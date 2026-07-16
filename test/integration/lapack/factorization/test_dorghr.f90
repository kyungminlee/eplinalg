! dorghr forms Q from a Hessenberg-reduction tau. Mirror dgehrd path.
program test_dorghr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgehrd, target_dorghr
    use ref_quad_lapack, only: dgehrd, dorghr
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorghr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A0, seed = 37001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), tau_ref(n-1), tau_got(n-1))
        A_ref = A0; A_got = A0
        call dgehrd(n, 1, n, A_ref, n, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A_ref, n, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgehrd(n, 1, n, A_got, n, tau_got, info)
        call dorghr(n, 1, n, A_ref, n, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dorghr(n, 1, n, A_ref, n, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dorghr(n, 1, n, A_got, n, tau_got, info)
        err = max_rel_err_mat(A_got, A_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dorghr
