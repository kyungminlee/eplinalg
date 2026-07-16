program test_blas_dsymv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dsymv_x, &
                              blas_upper, blas_lower
    use ref_quad_xblas, only: ref_blas_dsymv_x
    implicit none

    integer, parameter :: cases(3) = [10, 64, 200]
    integer :: c, n, uplo
    real(ep), allocatable :: a(:,:), x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('blas_dsymv_x', target_name)

    do c = 1, size(cases)
        n = cases(c)
        call gen_matrix_quad(n, n, a, seed = 100 + 13 * c)
        call gen_vector_quad(n, x,    seed = 200 + 13 * c)
        call gen_vector_quad(n, y_ref, seed = 300 + 13 * c)
        allocate(y_got(n))
        alpha = 1.7_ep
        beta  = -0.3_ep

        do uplo = 0, 1   ! 0=upper, 1=lower
            y_got = y_ref
            call ref_blas_dsymv_x(merge(blas_lower, blas_upper, uplo == 1), &
                                  n, alpha, a, n, x, 1, beta, y_ref, 1)
            call target_blas_dsymv_x(merge(blas_lower, blas_upper, uplo == 1), &
                                     n, alpha, a, n, x, 1, beta, y_got, 1)
            err = max_rel_err_vec(y_got, y_ref)
            tol = 16.0_ep * (2.0_ep * real(n, ep) * real(n, ep)) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ' uplo=', uplo
            call report_case(trim(label), err, tol)
            ! Restore y_ref baseline for next uplo iteration
            call gen_vector_quad(n, y_ref, seed = 300 + 13 * c)
        end do
        deallocate(a, x, y_ref, y_got)
    end do

    call report_finalize()
end program test_blas_dsymv_x
