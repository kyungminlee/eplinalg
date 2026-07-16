program test_blas_zhpmv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zhpmv_x, &
                              blas_upper
    use ref_quad_xblas, only: ref_blas_zhpmv_x
    implicit none
    integer, parameter :: n = 64
    integer, parameter :: lap = n * (n + 1) / 2
    complex(ep), allocatable :: ap(:), x(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    integer :: j, idx
    call report_init('blas_zhpmv_x', target_name)
    call gen_vector_complex(lap, ap, seed = 11)
    ! Diagonal must be real for hermitian (upper packed: a(j,j) at idx j+j*(j-1)/2)
    do j = 1, n
        idx = j + j * (j - 1) / 2
        ap(idx) = cmplx(real(ap(idx), ep), 0.0_ep, ep)
    end do
    call gen_vector_complex(n, x, seed = 12)
    call gen_vector_complex(n, y_ref, seed = 13)
    allocate(y_got(n)); y_got = y_ref
    alpha = (1.7_ep, -0.4_ep); beta = (-0.3_ep, 0.6_ep)
    call ref_blas_zhpmv_x(blas_upper, n, alpha, ap, x, 1, beta, y_ref, 1)
    call target_blas_zhpmv_x(blas_upper, n, alpha, ap, x, 1, beta, y_got, 1)
    err = max_rel_err_vec_z(y_got, y_ref)
    tol = 16.0_ep * (4.0_ep * real(n, ep) * real(n, ep)) * target_eps
    call report_case('n=64', err, tol)
    call report_finalize()
end program
