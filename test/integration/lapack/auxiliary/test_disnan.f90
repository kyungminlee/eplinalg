! disnan: NaN test for double precision.
program test_disnan
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use target_lapack,   only: target_name, target_eps, target_disnan
    use ref_quad_lapack, only: disnan
    implicit none

    real(ep) :: vals(7), zero
    logical :: ref_val, got_val
    integer :: i
    character(len=48) :: label

    zero = 0.0_ep
    vals(1) = 0.0_ep
    vals(2) = 1.0_ep
    vals(3) = -1.0_ep
    vals(4) = 1.0e30_ep
    vals(5) = huge(1.0_ep)
    vals(6) = tiny(1.0_ep)
    vals(7) = zero / zero

    call report_init('disnan', target_name)
    do i = 1, size(vals)
        ref_val = disnan(vals(i))
        got_val = target_disnan(vals(i))
        write(label, '(a,i0)') 'case=', i
        if (ref_val .eqv. got_val) then
            call report_case(trim(label), 0.0_ep, target_eps)
        else
            call report_case(trim(label), 1.0_ep, target_eps)
        end if
    end do
    call report_finalize()
end program test_disnan
