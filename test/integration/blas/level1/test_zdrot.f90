program test_zdrot
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zdrot
    use ref_quad_blas, only: zdrot
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x0(:), y0(:), x_ref(:), y_ref(:), x_got(:), y_got(:)
    real(ep) :: c, s, theta, err_x, err_y, tol
    character(len=32) :: label

    call report_init('zdrot', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x0, seed = 511 + 13 * i)
        call gen_vector_complex(n, y0, seed = 611 + 13 * i)
        theta = 0.3_ep + 0.1_ep * real(i, ep)
        c = cos(theta); s = sin(theta)
        allocate(x_ref(n), y_ref(n), x_got(n), y_got(n))
        x_ref = x0; y_ref = y0
        x_got = x0; y_got = y0
        call zdrot(n, x_ref, 1, y_ref, 1, c, s)
        call target_zdrot(n, x_got, 1, y_got, 1, c, s)
        err_x = max_rel_err_vec_z(x_got, x_ref)
        err_y = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), max(err_x, err_y), tol)
        deallocate(x_ref, y_ref, x_got, y_got)
    end do
    call report_finalize()
end program test_zdrot
