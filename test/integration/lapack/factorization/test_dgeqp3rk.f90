! dgeqp3rk: rank-revealing QR with column pivoting and truncation.
program test_dgeqp3rk
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgeqp3rk
    use ref_quad_lapack, only: dgeqp3rk
    implicit none

    integer, parameter :: ms(*) = [16, 24]
    integer, parameter :: ns(*) = [8,  12]
    integer, parameter :: nrhs  = 2
    integer :: i, m, n, kmax, info_ref, info_got, lwork, ncols
    integer :: K_ref, K_got
    real(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: tau_ref(:), tau_got(:), work(:)
    integer, allocatable :: jpiv_ref(:), jpiv_got(:), iwork(:)
    real(ep) :: wopt(1), err, tol, m2k_ref, m2k_got, rm2k_ref, rm2k_got
    character(len=48) :: label

    call report_init('dgeqp3rk', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kmax = n; ncols = n + nrhs
        call gen_matrix_quad(m, ncols, A0, seed = 23601 + 73 * i)
        allocate(A_ref(m, ncols), A_got(m, ncols), tau_ref(min(m, n)), tau_got(min(m, n)), &
                 jpiv_ref(n), jpiv_got(n), iwork(n))
        A_ref = A0; A_got = A0
        jpiv_ref = 0; jpiv_got = 0
        call dgeqp3rk(m, n, nrhs, kmax, -1.0_ep, -1.0_ep, A_ref, m, &
                      K_ref, m2k_ref, rm2k_ref, jpiv_ref, tau_ref, &
                      wopt, -1, iwork, info_ref)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqp3rk(m, n, nrhs, kmax, -1.0_ep, -1.0_ep, A_ref, m, &
                      K_ref, m2k_ref, rm2k_ref, jpiv_ref, tau_ref, &
                      work, lwork, iwork, info_ref)
        deallocate(work)
        call target_dgeqp3rk(m, n, nrhs, kmax, -1.0_ep, -1.0_ep, A_got, m, &
                             K_got, m2k_got, rm2k_got, jpiv_got, tau_got, info_got)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_vec(tau_got, tau_ref))
        if (K_ref /= K_got .or. any(jpiv_ref /= jpiv_got)) err = max(err, 1.0_ep)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, tau_ref, tau_got, jpiv_ref, jpiv_got, iwork)
    end do
    call report_finalize()
end program test_dgeqp3rk
