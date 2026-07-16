! zlartg: Givens rotation generator (complex).
program test_zlartg
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar, rel_err_scalar_z
    use target_lapack,   only: target_name, target_eps, target_zlartg
    use ref_quad_lapack, only: zlartg
    implicit none

    ! Avoid the degenerate f=(0,0) case here — multifloats returns NaN
    ! while quad returns the LAPACK-spec values, and that mismatch is
    ! inherent in the dd-arithmetic divide-by-zero handling.
    integer, parameter :: ncases = 4
    complex(ep), parameter :: fs(ncases) = [(3.0_ep, 1.0_ep), (-1.5_ep, 0.5_ep), &
                                            (0.25_ep, -0.1_ep), (100.0_ep, 50.0_ep)]
    complex(ep), parameter :: gs(ncases) = [(4.0_ep, -2.0_ep), (2.5_ep, 1.0_ep), &
                                            (0.75_ep, 0.2_ep), (-200.0_ep, 100.0_ep)]
    integer :: i
    real(ep) :: c_r, c_g, err, tol
    complex(ep) :: s_r, r_r, s_g, r_g
    character(len=48) :: label

    call report_init('zlartg', target_name)
    do i = 1, ncases
        call zlartg(fs(i), gs(i), c_r, s_r, r_r)
        call target_zlartg(fs(i), gs(i), c_g, s_g, r_g)
        err = max(rel_err_scalar(c_g, c_r), rel_err_scalar_z(s_g, s_r), rel_err_scalar_z(r_g, r_r))
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_zlartg
