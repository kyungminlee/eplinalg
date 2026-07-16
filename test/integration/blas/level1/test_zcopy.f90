program test_zcopy
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zcopy
    use ref_quad_blas, only: zcopy
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x(:), y_ref(:), y_got(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zcopy', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 121 + 17 * i)
        allocate(y_ref(n), y_got(n))
        y_ref = (0.0_ep, 0.0_ep); y_got = (0.0_ep, 0.0_ep)
        call zcopy(n, x, 1, y_ref, 1)
        call target_zcopy(n, x, 1, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zcopy
