program test_blas_zaxpby_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zaxpby_x
    use ref_quad_xblas, only: ref_blas_zaxpby_x
    implicit none
    integer, parameter :: cases(3) = [10, 100, 1000]
    integer :: i, n
    complex(ep), allocatable :: x(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=32) :: label
    call report_init('blas_zaxpby_x', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 100 + i)
        call gen_vector_complex(n, y_ref, seed = 200 + i)
        allocate(y_got(n)); y_got = y_ref
        alpha = (1.7_ep, -0.4_ep); beta = (-0.3_ep, 0.6_ep)
        call ref_blas_zaxpby_x(n, alpha, x, 1, beta, y_ref, 1)
        call target_blas_zaxpby_x(n, alpha, x, 1, beta, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 16.0_ep * 8.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x, y_ref, y_got)
    end do
    call report_finalize()
end program
