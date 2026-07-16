program test_zgelqf
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgelqf
    use ref_quad_lapack, only: zgelqf
    implicit none

    integer, parameter :: ms(*) = [8, 16, 32]
    integer, parameter :: ns(*) = [16, 24, 48]
    integer :: i, m, n, info, lwork, kmn
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err_a, err_t, tol
    character(len=48) :: label

    call report_init('zgelqf', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); kmn = min(m, n)
        call gen_matrix_complex(m, n, A0, seed = 40001 + 47 * i)
        allocate(A_ref(m,n), A_got(m,n), tau_ref(kmn), tau_got(kmn))
        A_ref = A0; A_got = A0
        call zgelqf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgelqf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgelqf(m, n, A_got, m, tau_got, info)
        err_a = max_rel_err_mat_z(A_got, A_ref)
        err_t = max_rel_err_vec_z(tau_got, tau_ref)
        tol   = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=A'
        call report_case(trim(label), err_a, tol)
        write(label, '(a,i0,a,i0,a)') 'm=', m, ',n=', n, ',out=tau'
        call report_case(trim(label), err_t, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zgelqf
