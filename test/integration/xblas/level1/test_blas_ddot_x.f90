program test_blas_ddot_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: rel_err_scalar
    use test_data,      only: gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_ddot_x, &
                              blas_no_conj => blas_no_conj
    use ref_quad_xblas, only: ref_blas_ddot_x
    implicit none

    integer, parameter :: cases(3) = [10, 100, 1000]
    integer :: i, n
    real(ep), allocatable :: x(:), y(:)
    real(ep) :: alpha, beta, ref_r, got_r, err, tol
    character(len=32) :: label

    call report_init('blas_ddot_x', target_name)

    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 42 + 7 * i)
        call gen_vector_quad(n, y, seed = 99 + 7 * i)
        alpha = 1.7_ep
        beta  = -0.3_ep

        ref_r = 0.5_ep
        got_r = 0.5_ep

        call ref_blas_ddot_x(blas_no_conj, n, alpha, x, 1, beta, y, 1, ref_r)
        call target_blas_ddot_x(blas_no_conj, n, alpha, x, 1, beta, y, 1, got_r)

        err = rel_err_scalar(got_r, ref_r)
        tol = 16.0_ep * (2.0_ep * real(n, ep)) * target_eps

        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
    end do

    call report_finalize()
end program test_blas_ddot_x
