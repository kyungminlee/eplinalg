program test_zdscal
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zdscal
    use ref_quad_blas, only: zdscal
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x0(:), x_ref(:), x_got(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('zdscal', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x0, seed = 421 + 13 * i)
        alpha = 0.75_ep + 0.05_ep * real(i, ep)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call zdscal(n, alpha, x_ref, 1)
        call target_zdscal(n, alpha, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_zdscal
