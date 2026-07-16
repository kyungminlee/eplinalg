! dlartgs: Givens rotation for SVD shifts.
program test_dlartgs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use target_lapack,   only: target_name, target_eps, target_dlartgs
    use ref_quad_lapack, only: dlartgs
    implicit none

    integer, parameter :: ncases = 5
    real(ep), parameter :: xs(ncases) = [3.0_ep, 1.5_ep, 5.0_ep,  10.0_ep, 0.5_ep]
    real(ep), parameter :: ys(ncases) = [4.0_ep, 2.5_ep, 1.0_ep, -2.0_ep,  0.7_ep]
    real(ep), parameter :: sg(ncases) = [1.0_ep, 0.8_ep, 2.0_ep, 0.5_ep,   0.1_ep]
    integer :: i
    real(ep) :: c_r, s_r, c_g, s_g, err, tol
    character(len=48) :: label

    call report_init('dlartgs', target_name)
    do i = 1, ncases
        call dlartgs(xs(i), ys(i), sg(i), c_r, s_r)
        call target_dlartgs(xs(i), ys(i), sg(i), c_g, s_g)
        err = max(rel_err_scalar(c_g, c_r), rel_err_scalar(s_g, s_r))
        tol = 64.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dlartgs
