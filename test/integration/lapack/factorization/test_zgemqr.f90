program test_zgemqr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgemqr
    use ref_quad_lapack, only: zgeqr, zgemqr
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [12, 20, 28]
    integer, parameter :: ncols = 5
    integer :: i, m, n, info, lwork, tsize
    complex(ep), allocatable :: A(:,:), T(:), work(:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: tquery(5), wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgemqr', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A, seed = 370061 + 47 * i)
        call zgeqr(m, n, A, m, tquery, -1, wopt, -1, info)
        tsize = max(1, int(real(tquery(1), ep)))
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(T(tsize), work(lwork))
        call zgeqr(m, n, A, m, T, tsize, work, lwork, info)
        deallocate(work)
        call gen_matrix_complex(m, ncols, C0, seed = 370071 + 47 * i)
        allocate(C_ref(m, ncols), C_got(m, ncols))
        C_ref = C0; C_got = C0
        call zgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_ref, m, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_ref, m, work, lwork, info)
        call target_zgemqr('L', 'N', m, ncols, n, A, m, T, tsize, C_got, m, info)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, T, C0, C_ref, C_got, work)
    end do
    call report_finalize()
end program test_zgemqr
