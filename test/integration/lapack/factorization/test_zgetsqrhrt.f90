! zgetsqrhrt: TSQR + Householder reconstruction (complex).
program test_zgetsqrhrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgetsqrhrt
    use ref_quad_lapack, only: zgetsqrhrt
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer :: i, m, n, mb1, nb1, nb2, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_ref(:,:), T_got(:,:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgetsqrhrt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb1 = 2*n; nb1 = 2; nb2 = 2  ! MB1 > N is required
        call gen_matrix_complex(m, n, A0, seed = 23551 + 71 * i)
        allocate(A_ref(m, n), A_got(m, n), T_ref(nb2, n), T_got(nb2, n))
        A_ref = A0; A_got = A0
        T_ref = (0.0_ep, 0.0_ep); T_got = (0.0_ep, 0.0_ep)
        call zgetsqrhrt(m, n, mb1, nb1, nb2, A_ref, m, T_ref, nb2, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgetsqrhrt(m, n, mb1, nb1, nb2, A_ref, m, T_ref, nb2, work, lwork, info)
        deallocate(work)
        call target_zgetsqrhrt(m, n, mb1, nb1, nb2, A_got, m, T_got, nb2, info)
        err = max(max_rel_err_mat_z(A_got, A_ref), max_rel_err_mat_z(T_got, T_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got)
    end do
    call report_finalize()
end program test_zgetsqrhrt
