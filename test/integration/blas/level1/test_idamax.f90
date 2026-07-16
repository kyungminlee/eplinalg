program test_idamax
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_idamax
    use ref_quad_blas, only: idamax
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n, ref, got
    real(ep), allocatable :: x(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('idamax', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 71 + 13 * i)
        ref = idamax(n, x, 1)
        got = target_idamax(n, x, 1)
        ! Index match: error is 0 on match, 1 otherwise. Tolerance 0.
        if (got == ref) then
            err = 0.0_ep
        else
            err = 1.0_ep
        end if
        tol = 0.0_ep
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        if (allocated(x)) deallocate(x)
    end do

    ! Tie-breaking: BLAS specifies the *first* index of the maximum.
    ! Random vectors hit ties with vanishing probability, so construct
    ! a plateau explicitly: x = [-1, 2, 3, 3, 3, -3, 1, 0, 3, 2].
    ! All four ``3`` / ``-3`` entries have the same |.|; idamax must
    ! return 3 (1-based) — the first.
    allocate(x(10))
    x = [-1.0_ep, 2.0_ep, 3.0_ep, 3.0_ep, 3.0_ep, &
         -3.0_ep, 1.0_ep, 0.0_ep, 3.0_ep, 2.0_ep]
    ref = idamax(10, x, 1)
    got = target_idamax(10, x, 1)
    if (got == ref) then
        err = 0.0_ep
    else
        err = 1.0_ep
    end if
    call report_case('plateau-tie', err, 0.0_ep)
    deallocate(x)

    call report_finalize()
end program test_idamax
