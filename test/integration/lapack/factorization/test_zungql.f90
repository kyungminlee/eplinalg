program test_zungql
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zgeqlf, target_zungql
    use ref_quad_lapack, only: zgeqlf, zungql
    implicit none

    integer, parameter :: ms(*) = [16, 48, 96]
    integer, parameter :: ns(*) = [8,  24, 64]
    integer :: i, m, n, info, lwork
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: tau_ref(:), tau_got(:), work(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zungql', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 64001 + 61 * i)
        allocate(A_ref(m,n), A_got(m,n), tau_ref(min(m,n)), tau_got(min(m,n)))
        A_ref = A0; A_got = A0
        call zgeqlf(m, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgeqlf(m, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zgeqlf(m, n, A_got, m, tau_got, info)
        call zungql(m, n, n, A_ref, m, tau_ref, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zungql(m, n, n, A_ref, m, tau_ref, work, lwork, info)
        deallocate(work)
        call target_zungql(m, n, n, A_got, m, tau_got, info)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, tau_ref, tau_got)
    end do
    call report_finalize()
end program test_zungql
