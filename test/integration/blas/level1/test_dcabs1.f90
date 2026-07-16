program test_dcabs1
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use target_blas,   only: target_name, target_eps, target_dcabs1
    use ref_quad_blas, only: dcabs1
    implicit none

    integer, parameter :: cases = 8
    integer :: i
    complex(ep) :: z
    real(ep) :: re, im, ref, got, err, tol
    character(len=32) :: label

    call report_init('dcabs1', target_name)
    do i = 1, cases
        re = real(i, ep) - 4.5_ep
        im = real(2 * i, ep) - 7.0_ep
        z = cmplx(re, im, ep)
        ref = dcabs1(z)
        got = target_dcabs1(z)
        err = rel_err_scalar(got, ref)
        tol = 4.0_ep * target_eps
        write(label, '(a,i0)') 'i=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dcabs1
