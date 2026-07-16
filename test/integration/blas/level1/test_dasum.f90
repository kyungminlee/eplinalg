program test_dasum
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dasum
    use ref_quad_blas, only: dasum
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call report_init('dasum', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 11 + 13 * i)
        ref = dasum(n, x, 1)
        got = target_dasum(n, x, 1)
        err = rel_err_scalar(got, ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dasum
