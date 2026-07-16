program test_blas_zsum_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: rel_err_scalar_z
    use test_data,      only: gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zsum_x
    use ref_quad_xblas, only: ref_blas_zsum_x
    implicit none
    integer, parameter :: cases(3) = [10, 100, 1000]
    integer :: i, n
    complex(ep), allocatable :: x(:)
    complex(ep) :: ref_s, got_s
    real(ep) :: err, tol
    character(len=32) :: label
    call report_init('blas_zsum_x', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 100 + i)
        call ref_blas_zsum_x(n, x, 1, ref_s)
        call target_blas_zsum_x(n, x, 1, got_s)
        err = rel_err_scalar_z(got_s, ref_s)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(x)
    end do
    call report_finalize()
end program
