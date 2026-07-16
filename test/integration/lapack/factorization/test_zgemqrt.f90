program test_zgemqrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgemqrt
    use ref_quad_lapack, only: zgeqrt, zgemqrt
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [12, 20, 28]
    integer, parameter :: ncols = 5
    integer, parameter :: nb    = 4
    integer :: i, m, n, info
    complex(ep), allocatable :: V(:,:), T(:,:), C0(:,:), C_ref(:,:), C_got(:,:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgemqrt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, V, seed = 380061 + 47 * i)
        allocate(T(nb, min(m, n)), work(nb*max(m, n)))
        call zgeqrt(m, n, nb, V, m, T, nb, work, info)
        deallocate(work)
        call gen_matrix_complex(m, ncols, C0, seed = 380071 + 47 * i)
        allocate(C_ref(m, ncols), C_got(m, ncols), work(ncols*nb))
        C_ref = C0; C_got = C0
        call zgemqrt('L', 'N', m, ncols, n, nb, V, m, T, nb, C_ref, m, work, info)
        call target_zgemqrt('L', 'N', m, ncols, n, nb, V, m, T, nb, C_got, m, info)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(V, T, C0, C_ref, C_got, work)
    end do
    call report_finalize()
end program test_zgemqrt
