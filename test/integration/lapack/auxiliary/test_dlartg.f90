! dlartg: Givens rotation generator.
program test_dlartg
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use target_lapack,   only: target_name, target_eps, target_dlartg
    use ref_quad_lapack, only: dlartg
    implicit none

    integer, parameter :: ncases = 5
    real(ep), parameter :: fs(ncases) = [3.0_ep, -1.5_ep, 0.25_ep,  100.0_ep,  0.0_ep]
    real(ep), parameter :: gs(ncases) = [4.0_ep,  2.5_ep, 0.75_ep, -200.0_ep, 7.0_ep]
    integer :: i
    real(ep) :: c_r, s_r, r_r, c_g, s_g, r_g, err, tol
    character(len=48) :: label

    call report_init('dlartg', target_name)
    do i = 1, ncases
        call dlartg(fs(i), gs(i), c_r, s_r, r_r)
        call target_dlartg(fs(i), gs(i), c_g, s_g, r_g)
        err = max(rel_err_scalar(c_g, c_r), rel_err_scalar(s_g, s_r), rel_err_scalar(r_g, r_r))
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dlartg
