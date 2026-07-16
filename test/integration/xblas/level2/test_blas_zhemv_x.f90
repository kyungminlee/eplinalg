program test_blas_zhemv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_matrix_complex, gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zhemv_x, &
                              blas_upper, blas_lower
    use ref_quad_xblas, only: ref_blas_zhemv_x
    implicit none

    integer, parameter :: cases(3) = [10, 64, 200]
    integer :: c, n, uplo, i
    complex(ep), allocatable :: a(:,:), x(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('blas_zhemv_x', target_name)

    do c = 1, size(cases)
        n = cases(c)
        call gen_matrix_complex(n, n, a, seed = 100 + 13 * c)
        ! Force diagonal to be real — Hermitian matrices have real diagonals.
        do i = 1, n
            a(i, i) = cmplx(real(a(i, i), ep), 0.0_ep, ep)
        end do
        call gen_vector_complex(n, x,    seed = 200 + 13 * c)
        call gen_vector_complex(n, y_ref, seed = 300 + 13 * c)
        allocate(y_got(n))
        alpha = (1.7_ep, -0.4_ep)
        beta  = (-0.3_ep, 0.6_ep)

        do uplo = 0, 1
            y_got = y_ref
            call ref_blas_zhemv_x(merge(blas_lower, blas_upper, uplo == 1), &
                                  n, alpha, a, n, x, 1, beta, y_ref, 1)
            call target_blas_zhemv_x(merge(blas_lower, blas_upper, uplo == 1), &
                                     n, alpha, a, n, x, 1, beta, y_got, 1)
            err = max_rel_err_vec_z(y_got, y_ref)
            tol = 16.0_ep * (4.0_ep * real(n, ep) * real(n, ep)) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ' uplo=', uplo
            call report_case(trim(label), err, tol)
            call gen_vector_complex(n, y_ref, seed = 300 + 13 * c)
        end do
        deallocate(a, x, y_ref, y_got)
    end do

    call report_finalize()
end program test_blas_zhemv_x
