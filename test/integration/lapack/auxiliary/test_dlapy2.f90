! dlapy2: hypotenuse sqrt(x^2+y^2) without spurious overflow.
program test_dlapy2
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use target_lapack,   only: target_name, target_eps, target_dlapy2
    use ref_quad_lapack, only: dlapy2
    implicit none

    integer, parameter :: ncases = 5
    real(ep), parameter :: xs(ncases) = [3.0_ep, -1.5_ep, 0.25_ep, 100.0_ep, 1.0e30_ep]
    real(ep), parameter :: ys(ncases) = [4.0_ep,  2.5_ep, 0.75_ep, -200.0_ep, 1.0e30_ep]
    integer :: i
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('dlapy2', target_name)
    do i = 1, ncases
        ref_val = dlapy2(xs(i), ys(i))
        got_val = target_dlapy2(xs(i), ys(i))
        err = rel_err_scalar(got_val, ref_val)
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dlapy2
