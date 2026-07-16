program test_zgehrd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgehrd
    use ref_quad_lapack, only: zgehrd
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), tau_ref(:), tau_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err_a, err_t, tol
    character(len=48) :: label

    call report_init('zgehrd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 46001 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), tau_ref(n-1), tau_got(n-1))
        A_ref = A0; A_got = A0
        call zgehrd(n, 1, n, A_ref, n, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgehrd(n, 1, n, A_ref, n, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgehrd(n, 1, n, A_got, n, tau_got, info)
        err_a = max_rel_err_mat_z(A_got, A_ref)
        err_t = max_rel_err_vec_z(tau_got, tau_ref)
        tol   = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a)') 'n=', n, ',out=A'
        call report_case(trim(label), err_a, tol)
        write(label, '(a,i0,a)') 'n=', n, ',out=tau'
        call report_case(trim(label), err_t, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zgehrd
