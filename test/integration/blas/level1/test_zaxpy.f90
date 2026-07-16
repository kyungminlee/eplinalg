program test_zaxpy
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zaxpy
    use ref_quad_blas, only: zaxpy
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zaxpy', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x,  seed = 211 + 13 * i)
        call gen_vector_complex(n, y0, seed = 311 + 13 * i)
        alpha = cmplx(0.7_ep + 0.1_ep * i, 0.3_ep - 0.05_ep * i, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0
        y_got = y0
        call zaxpy(n, alpha, x, 1, y_ref, 1)
        call target_zaxpy(n, alpha, x, 1, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * target_eps  ! complex axpy ≈ 8 FLOPs per element
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zaxpy
