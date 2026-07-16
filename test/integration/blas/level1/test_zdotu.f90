program test_zdotu
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zdotu
    use ref_quad_blas, only: zdotu
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x(:), y(:)
    complex(ep) :: ref, got
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('zdotu', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 231 + 13 * i)
        call gen_vector_complex(n, y, seed = 331 + 13 * i)
        ref = zdotu(n, x, 1, y, 1)
        got = target_zdotu(n, x, 1, y, 1)
        err = rel_err_scalar_z(got, ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_zdotu
