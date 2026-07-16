! dgesvdx: SVD with index/value selection (range='A', 'V', or 'I').
program test_dgesvdx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgesvdx
    use ref_quad_lapack, only: dgesvdx
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [8,  16]
    integer :: i, m, n, info, ns_r, ns_g, lwork
    real(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:), S_r(:), S_g(:), U(:,:), Vt(:,:), work(:)
    integer, allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgesvdx', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 25101 + 73 * i)
        allocate(A_r(m, n), A_g(m, n), S_r(min(m, n)), S_g(min(m, n)), &
                 U(m, min(m, n)), Vt(min(m, n), n), iwork(12*min(m, n)))
        A_r = A0; A_g = A0
        call dgesvdx('N', 'N', 'A', m, n, A_r, m, 0.0_ep, 0.0_ep, 0, 0, &
                     ns_r, S_r, U, m, Vt, min(m, n), wopt, -1, iwork, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgesvdx('N', 'N', 'A', m, n, A_r, m, 0.0_ep, 0.0_ep, 0, 0, &
                     ns_r, S_r, U, m, Vt, min(m, n), work, lwork, iwork, info)
        deallocate(work)
        call target_dgesvdx('N', 'N', 'A', m, n, A_g, m, 0.0_ep, 0.0_ep, 0, 0, &
                            ns_g, S_g, U, m, Vt, min(m, n), info)
        err = max_rel_err_vec(S_g(1:ns_g), S_r(1:ns_r))
        if (ns_r /= ns_g) err = max(err, 1.0_ep)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, S_r, S_g, U, Vt, iwork)
    end do
    call report_finalize()
end program test_dgesvdx
