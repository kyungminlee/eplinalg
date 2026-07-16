program test_dcopy
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dcopy
    use ref_quad_blas, only: dcopy
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x(:), y_ref(:), y_got(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('dcopy', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 31 + 13 * i)
        allocate(y_ref(n), y_got(n))
        y_ref = -1.0_ep
        y_got = -1.0_ep
        call dcopy(n, x, 1, y_ref, 1)
        call target_dcopy(n, x, 1, y_got, 1)
        err = max_rel_err_vec(y_got, y_ref)
        tol = target_eps  ! copy is exact
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_dcopy
