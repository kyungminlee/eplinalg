program test_zrotg
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar, rel_err_scalar_z
    use target_blas,   only: target_name, target_eps, target_zrotg
    use ref_quad_blas, only: zrotg
    implicit none

    integer, parameter :: cases = 6
    integer :: i
    complex(ep) :: a0, b, a_ref, a_got, s_ref, s_got
    real(ep) :: c_ref, c_got, err, tol
    character(len=32) :: label

    call report_init('zrotg', target_name)
    do i = 1, cases
        a0 = cmplx(real(i, ep) * 0.5_ep + 0.1_ep, real(i, ep) * 0.3_ep - 0.2_ep, ep)
        b  = cmplx(0.7_ep - 0.05_ep * real(i, ep), 0.4_ep + 0.1_ep * real(i, ep), ep)
        a_ref = a0; a_got = a0
        call zrotg(a_ref, b, c_ref, s_ref)
        call target_zrotg(a_got, b, c_got, s_got)
        err = max(rel_err_scalar_z(a_got, a_ref), &
                  rel_err_scalar(c_got, c_ref),  &
                  rel_err_scalar_z(s_got, s_ref))
        tol = 32.0_ep * target_eps
        write(label, '(a,i0)') 'i=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_zrotg
