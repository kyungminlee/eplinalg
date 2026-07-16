! dlapy3: 3-arg hypot sqrt(x^2+y^2+z^2).
program test_dlapy3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use target_lapack,   only: target_name, target_eps, target_dlapy3
    use ref_quad_lapack, only: dlapy3
    implicit none

    integer, parameter :: ncases = 5
    real(ep), parameter :: xs(ncases) = [1.0_ep, -2.0_ep, 0.1_ep, 100.0_ep, 1.0e15_ep]
    real(ep), parameter :: ys(ncases) = [2.0_ep,  3.0_ep, 0.2_ep, 200.0_ep, 1.0e15_ep]
    real(ep), parameter :: zs(ncases) = [3.0_ep, -4.0_ep, 0.3_ep, 300.0_ep, 1.0e15_ep]
    integer :: i
    real(ep) :: ref_val, got_val, err, tol
    character(len=48) :: label

    call report_init('dlapy3', target_name)
    do i = 1, ncases
        ref_val = dlapy3(xs(i), ys(i), zs(i))
        got_val = target_dlapy3(xs(i), ys(i), zs(i))
        err = rel_err_scalar(got_val, ref_val)
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dlapy3
