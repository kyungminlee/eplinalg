program test_daxpy
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_daxpy
    use ref_quad_blas, only: daxpy
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x(:), y0(:), y_ref(:), y_got(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('daxpy', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x,  seed = 21 + 13 * i)
        call gen_vector_quad(n, y0, seed = 121 + 13 * i)
        alpha = real(0.7_ep + 0.1_ep * i, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0
        y_got = y0
        call daxpy(n, alpha, x, 1, y_ref, 1)
        call target_daxpy(n, alpha, x, 1, y_got, 1)
        err = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * 2.0_ep * target_eps  ! axpy is 2 FLOPs per element
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_daxpy
