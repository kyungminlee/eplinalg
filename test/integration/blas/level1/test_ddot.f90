program test_ddot
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: rel_err_scalar
    use test_data,     only: gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_ddot
    use ref_quad_blas, only: ddot
    implicit none

    integer, parameter :: cases(*) = [10, 100, 1000, 10000]
    integer :: i, n
    real(ep), allocatable :: x(:), y(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call report_init('ddot', target_name)

    do i = 1, size(cases)
        n = cases(i)

        call gen_vector_quad(n, x, seed = 42 + 7 * i)
        call gen_vector_quad(n, y, seed = 99 + 7 * i)

        ref = ddot(n, x, 1, y, 1)
        got = target_ddot(n, x, 1, y, 1)

        err = rel_err_scalar(got, ref)
        ! Tolerance: 2*n FLOPs * eps * safety factor 16. The eps used
        ! is the *target's* epsilon — for kind16, eps_quad ≈ 2e-34, so
        ! a million-element dot product still fits comfortably.
        tol = 16.0_ep * (2.0_ep * real(n, ep)) * target_eps

        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
    end do

    call report_finalize()
end program test_ddot
