! dladiv: complex division (a+bi)/(c+di) → (p+qi).
program test_dladiv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar
    use target_lapack,   only: target_name, target_eps, target_dladiv
    use ref_quad_lapack, only: dladiv
    implicit none

    integer, parameter :: ncases = 4
    real(ep), parameter :: as(ncases) = [3.0_ep, -1.5_ep, 1.0_ep, 100.0_ep]
    real(ep), parameter :: bs(ncases) = [2.0_ep,  1.0_ep, 0.0_ep, 200.0_ep]
    real(ep), parameter :: cs(ncases) = [1.0_ep,  0.5_ep, 2.0_ep, 50.0_ep]
    real(ep), parameter :: ds(ncases) = [4.0_ep, -2.0_ep, 1.0_ep, 25.0_ep]
    integer :: i
    real(ep) :: p_r, q_r, p_g, q_g, err, tol
    character(len=48) :: label

    call report_init('dladiv', target_name)
    do i = 1, ncases
        call dladiv(as(i), bs(i), cs(i), ds(i), p_r, q_r)
        call target_dladiv(as(i), bs(i), cs(i), ds(i), p_g, q_g)
        err = max(rel_err_scalar(p_g, p_r), rel_err_scalar(q_g, q_r))
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_dladiv
