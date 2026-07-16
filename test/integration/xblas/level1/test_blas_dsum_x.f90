program test_blas_dsum_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: rel_err_scalar
    use test_data,      only: gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dsum_x
    use ref_quad_xblas, only: ref_blas_dsum_x
    implicit none
    integer, parameter :: cases(3) = [10, 100, 1000]
    integer :: i, n
    real(ep), allocatable :: x(:)
    real(ep) :: ref_s, got_s, err, tol
    character(len=32) :: label
    call report_init('blas_dsum_x', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_quad(n, x, seed = 100 + i)
        call ref_blas_dsum_x(n, x, 1, ref_s)
        call target_blas_dsum_x(n, x, 1, got_s)
        err = rel_err_scalar(got_s, ref_s)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x)
    end do
    call report_finalize()
end program
