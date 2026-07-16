program test_zswap
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zswap
    use ref_quad_blas, only: zswap
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x0(:), y0(:), x_ref(:), y_ref(:), x_got(:), y_got(:)
    real(ep) :: err_x, err_y, tol
    character(len=32) :: label

    call report_init('zswap', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x0, seed = 221 + 11 * i)
        call gen_vector_complex(n, y0, seed = 321 + 11 * i)
        allocate(x_ref(n), y_ref(n), x_got(n), y_got(n))
        x_ref = x0; y_ref = y0
        x_got = x0; y_got = y0
        call zswap(n, x_ref, 1, y_ref, 1)
        call target_zswap(n, x_got, 1, y_got, 1)
        err_x = max_rel_err_vec_z(x_got, x_ref)
        err_y = max_rel_err_vec_z(y_got, y_ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), max(err_x, err_y), tol)
        deallocate(x_ref, y_ref, x_got, y_got)
    end do
    call report_finalize()
end program test_zswap
