program test_dasqrtb
    use prec_kinds,           only: ep
    use compare,              only: rel_err_scalar
    use ptzblas_prec_report,  only: report_init, report_case, report_finalize
    use ptzblas_ref_quad_aux, only: ref_dasqrtb
    use target_ptzblas,       only: target_name, target_eps, target_dasqrtb
    implicit none

    integer, parameter :: ncases = 6
    real(ep), parameter :: a_set(ncases) = [1.0_ep, -2.5_ep, 0.7_ep, 100.0_ep, -0.001_ep, 3.5_ep]
    real(ep), parameter :: b_set(ncases) = [4.0_ep,  2.0_ep, 0.25_ep, 1.0e-6_ep, 1.0e6_ep,  0.0_ep]
    integer :: i
    real(ep) :: c_got, c_ref, err, tol
    character(len=32) :: label

    call report_init('dasqrtb', target_name, 0)
    do i = 1, ncases
        call target_dasqrtb(a_set(i), b_set(i), c_got)
        call ref_dasqrtb(a_set(i), b_set(i), c_ref)
        err = rel_err_scalar(c_got, c_ref)
        tol = 8.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dasqrtb
