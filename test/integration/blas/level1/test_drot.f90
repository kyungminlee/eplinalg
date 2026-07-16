program test_drot
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_drot
    use ref_quad_blas, only: drot
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x0(:), y0(:)
    real(ep), allocatable :: x_ref(:), y_ref(:), x_got(:), y_got(:)
    real(ep) :: c, s, err_x, err_y, tol, theta
    character(len=32) :: label

    call report_init('drot', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x0, seed = 81 + 13 * i)
        call gen_vector_quad(n, y0, seed = 181 + 13 * i)
        ! Use a Givens rotation with a non-trivial angle.
        theta = 0.3_ep + 0.1_ep * i
        c = cos(theta)
        s = sin(theta)
        allocate(x_ref(n), y_ref(n), x_got(n), y_got(n))
        x_ref = x0; y_ref = y0
        x_got = x0; y_got = y0
        call drot(n, x_ref, 1, y_ref, 1, c, s)
        call target_drot(n, x_got, 1, y_got, 1, c, s)
        err_x = max_rel_err_vec(x_got, x_ref)
        err_y = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * 4.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), max(err_x, err_y), tol)
        deallocate(x_ref, y_ref, x_got, y_got)
    end do
    call report_finalize()
end program test_drot
