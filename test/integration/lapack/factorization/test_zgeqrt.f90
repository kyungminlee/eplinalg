! zgeqrt: recursive blocked QR with explicit T. Compares the full
! overwritten A and the upper triangle of each NB-by-NB T block; the
! strict-lower of T is workspace per LAPACK convention and is not
! compared.
program test_zgeqrt
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgeqrt
    use ref_quad_lapack, only: zgeqrt
    implicit none

    integer, parameter :: ms(*) = [16, 32, 48]
    integer, parameter :: ns(*) = [12, 20, 28]
    integer, parameter :: nb    = 4
    integer :: i, m, n, info, j, kmin, ii, ib, kk
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), T_ref(:,:), T_got(:,:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgeqrt', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 380011 + 47 * i)
        allocate(A_ref(m, n), A_got(m, n), T_ref(nb, min(m, n)), T_got(nb, min(m, n)), work(nb*n))
        A_ref = A0; A_got = A0
        call zgeqrt(m, n, nb, A_ref, m, T_ref, nb, work, info)
        call target_zgeqrt(m, n, nb, A_got, m, T_got, nb, info)
        kmin = min(m, n)
        err = max_rel_err_mat_z(A_got, A_ref)
        do ii = 1, kmin, nb
            ib = min(nb, kmin - ii + 1)
            do kk = 0, ib - 1
                do j = 1, kk + 1
                    err = max(err, abs(T_got(j, ii+kk) - T_ref(j, ii+kk)) &
                                   / max(abs(T_ref(j, ii+kk)), 1.0e-50_ep))
                end do
            end do
        end do
        tol = 16.0_ep * real(max(m, n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_ref, A_got, T_ref, T_got, work)
    end do
    call report_finalize()
end program test_zgeqrt
