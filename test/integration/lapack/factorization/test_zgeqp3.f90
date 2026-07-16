program test_zgeqp3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeqp3
    use ref_quad_lapack, only: zgeqp3
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [8,  24, 32]
    integer :: i, m, n, info, lwork, kmn, j, k
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    real(ep),    allocatable :: rwork(:)
    integer,     allocatable :: jpvt_ref(:), jpvt_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgeqp3', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kmn = min(m, n)
        call gen_matrix_complex(m, n, A0, seed = 49001 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n), tau_ref(kmn), tau_got(kmn))
        allocate(jpvt_ref(n), jpvt_got(n), rwork(2*n))
        A_ref = A0; A_got = A0
        jpvt_ref = 0; jpvt_got = 0
        call zgeqp3(m, n, A_ref, m, jpvt_ref, tau_ref, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgeqp3(m, n, A_ref, m, jpvt_ref, tau_ref, work, lwork, rwork, info)
        deallocate(work)
        call target_zgeqp3(m, n, A_got, m, jpvt_got, tau_got, info)
        do j = 1, n
            do k = min(j+1, m), m
                A_ref(k, j) = (0.0_ep, 0.0_ep); A_got(k, j) = (0.0_ep, 0.0_ep)
            end do
        end do
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got, jpvt_ref, jpvt_got, rwork)
    end do
    call report_finalize()
end program test_zgeqp3
