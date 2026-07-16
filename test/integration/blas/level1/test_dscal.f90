program test_dscal
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dscal
    use ref_quad_blas, only: dscal
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x0(:), x_ref(:), x_got(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call report_init('dscal', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x0, seed = 51 + 13 * i)
        alpha = real(0.5_ep + 0.13_ep * i, ep)
        allocate(x_ref(n), x_got(n))
        x_ref = x0
        x_got = x0
        call dscal(n, alpha, x_ref, 1)
        call target_dscal(n, alpha, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_dscal
