! dgeqrfp: blocked QR with non-negative R diagonal.
program test_dgeqrfp
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgeqrfp
    use ref_quad_lapack, only: dgeqrfp
    implicit none

    integer, parameter :: ms(*) = [16, 48]
    integer, parameter :: ns(*) = [8,  24]
    integer :: i, m, n, info, lwork
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgeqrfp', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 23101 + 73 * i)
        allocate(A_ref(m, n), A_got(m, n), tau_ref(n), tau_got(n))
        A_ref = A0; A_got = A0
        call dgeqrfp(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrfp(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_dgeqrfp(m, n, A_got, m, tau_got, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_vec(tau_got, tau_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_dgeqrfp
