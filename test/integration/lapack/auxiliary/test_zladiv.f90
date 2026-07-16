! zladiv: complex division x/y returning complex.
program test_zladiv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: rel_err_scalar_z
    use target_lapack,   only: target_name, target_eps, target_zladiv
    use ref_quad_lapack, only: zladiv
    implicit none

    integer, parameter :: ncases = 4
    complex(ep), parameter :: xs(ncases) = [(3.0_ep, 2.0_ep), (-1.5_ep, 1.0_ep), &
                                            (1.0_ep, 0.0_ep), (100.0_ep, 200.0_ep)]
    complex(ep), parameter :: ys(ncases) = [(1.0_ep, 4.0_ep), (0.5_ep, -2.0_ep), &
                                            (2.0_ep, 1.0_ep), (50.0_ep, 25.0_ep)]
    integer :: i
    complex(ep) :: ref_val, got_val
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zladiv', target_name)
    do i = 1, ncases
        ref_val = zladiv(xs(i), ys(i))
        got_val = target_zladiv(xs(i), ys(i))
        err = rel_err_scalar_z(got_val, ref_val)
        tol = 16.0_ep * target_eps
        write(label, '(a,i0)') 'case=', i
        call report_case(trim(label), err, tol)
    end do
    call report_finalize()
end program test_zladiv
