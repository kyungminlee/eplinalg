program test_dvvdot
    use prec_kinds,            only: ep
    use compare,               only: rel_err_scalar
    use ptzblas_prec_report,   only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_blas, only: ddot
    use test_data,             only: gen_vector_quad
    use target_ptzblas,        only: target_name, target_eps, target_dvvdot
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer :: i, n
    real(ep), allocatable :: x(:), y(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call report_init('dvvdot', target_name, 0)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 700 + 11 * i)
        call gen_vector_quad(n, y, seed = 800 + 13 * i)
        got = 0.0_ep
        call target_dvvdot(n, got, x, 1, y, 1)
        ref = ddot(n, x, 1, y, 1)
        err = rel_err_scalar(got, ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x, y)
    end do
    call report_finalize()
end program test_dvvdot
