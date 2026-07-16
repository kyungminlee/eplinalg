program test_drotg
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use target_blas,   only: target_name, target_eps, target_drotg
    use ref_quad_blas, only: drotg
    implicit none

    integer, parameter :: ncases = 5
    integer :: i
    real(ep) :: a_ref, b_ref, c_ref, s_ref
    real(ep) :: a_got, b_got, c_got, s_got
    real(ep) :: a0, b0, err, tol
    character(len=32) :: label

    call report_init('drotg', target_name)
    do i = 1, ncases
        ! Vary a/b magnitudes to exercise different branches.
        a0 = real(i, ep) * 1.5_ep
        b0 = real(ncases - i + 1, ep) * 0.7_ep
        a_ref = a0; b_ref = b0
        a_got = a0; b_got = b0
        call drotg(a_ref, b_ref, c_ref, s_ref)
        call target_drotg(a_got, b_got, c_got, s_got)
        err = max(rel_err_scalar(a_got, a_ref), &
                  rel_err_scalar(b_got, b_ref), &
                  rel_err_scalar(c_got, c_ref), &
                  rel_err_scalar(s_got, s_ref))
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_drotg
