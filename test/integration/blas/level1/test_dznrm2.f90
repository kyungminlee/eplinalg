program test_dznrm2
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_dznrm2
    use ref_quad_blas, only: dznrm2
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    complex(ep), allocatable :: x(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call report_init('dznrm2', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 71 + 13 * i)
        ref = dznrm2(n, x, 1)
        got = target_dznrm2(n, x, 1)
        err = rel_err_scalar(got, ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dznrm2
