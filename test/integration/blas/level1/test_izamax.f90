program test_izamax
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_izamax
    use ref_quad_blas, only: izamax
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n, ref, got
    complex(ep), allocatable :: x(:)
    real(ep) :: err
    character(len=32) :: label

    call report_init('izamax', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 91 + 13 * i)
        ref = izamax(n, x, 1)
        got = target_izamax(n, x, 1)
        if (got == ref) then
            err = 0.0_ep
        else
            err = 1.0_ep
        end if
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, 0.0_ep)
        if (allocated(x)) deallocate(x)
    end do

    ! Tie-breaking: izamax uses |re|+|im| (cabs1) and returns the FIRST
    ! index of the maximum. Build a plateau where indices 3, 5, 8 all
    ! tie at |re|+|im| = 5 — the result must be 3.
    allocate(x(10))
    x = [cmplx(1.0_ep,  0.0_ep, ep), cmplx( 0.0_ep, 2.0_ep, ep), &
         cmplx(3.0_ep,  2.0_ep, ep), cmplx( 1.0_ep, 1.0_ep, ep), &
         cmplx(4.0_ep, -1.0_ep, ep), cmplx(-2.0_ep, 2.0_ep, ep), &
         cmplx(0.0_ep,  0.0_ep, ep), cmplx(-3.0_ep, 2.0_ep, ep), &
         cmplx(1.0_ep,  3.0_ep, ep), cmplx( 2.0_ep, 1.0_ep, ep)]
    ref = izamax(10, x, 1)
    got = target_izamax(10, x, 1)
    if (got == ref) then
        err = 0.0_ep
    else
        err = 1.0_ep
    end if
    call report_case('plateau-tie', err, 0.0_ep)
    deallocate(x)

    call report_finalize()
end program test_izamax
