! zgeqr2p: unblocked QR with non-negative R diagonal (complex).
program test_zgeqr2p
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeqr2p
    use ref_quad_lapack, only: zgeqr2p
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [8,  16]
    integer :: i, m, n, info
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgeqr2p', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 23051 + 71 * i)
        allocate(A_ref(m, n), A_got(m, n), tau_ref(n), tau_got(n), work(n))
        A_ref = A0; A_got = A0
        call zgeqr2p(m, n, A_ref, m, tau_ref, work, info)
        call target_zgeqr2p(m, n, A_got, m, tau_got, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_vec_z(tau_got, tau_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got, work)
    end do
    call report_finalize()
end program test_zgeqr2p
