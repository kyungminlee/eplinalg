program test_zgerc
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zgerc
    use ref_quad_blas, only: zgerc
    implicit none

    integer, parameter :: ms(*) = [5, 50, 100]
    integer, parameter :: ns(*) = [7, 60, 120]
    integer :: i, m, n
    complex(ep), allocatable :: A0(:,:), x(:), y(:)
    complex(ep), allocatable :: A_ref(:,:), A_got(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zgerc', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_complex(m, n, A0, seed = 771 + 19 * i)
        call gen_vector_complex(m, x,   seed = 781 + 19 * i)
        call gen_vector_complex(n, y,   seed = 791 + 19 * i)
        alpha = cmplx(0.4_ep, 0.15_ep, ep)
        allocate(A_ref(m, n), A_got(m, n))
        A_ref = A0
        A_got = A0
        call zgerc(m, n, alpha, x, 1, y, 1, A_ref, m)
        call target_zgerc(m, n, alpha, x, 1, y, 1, A_got, m)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 32.0_ep * 8.0_ep * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got)
    end do
    call report_finalize()
end program test_zgerc
