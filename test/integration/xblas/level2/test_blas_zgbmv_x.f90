program test_blas_zgbmv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec_z
    use test_data,      only: gen_matrix_complex, gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zgbmv_x, &
                              blas_no_trans
    use ref_quad_xblas, only: ref_blas_zgbmv_x
    implicit none

    integer, parameter :: m = 64, n = 64, kl = 3, ku = 2
    integer, parameter :: lda = kl + ku + 1
    complex(ep), allocatable :: a(:,:), x(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('blas_zgbmv_x', target_name)

    call gen_matrix_complex(lda, n, a, seed = 1011)
    call gen_vector_complex(n, x,    seed = 1012)
    call gen_vector_complex(m, y_ref, seed = 1013)
    allocate(y_got(m))
    y_got = y_ref
    alpha = (1.7_ep, -0.4_ep)
    beta  = (-0.3_ep, 0.6_ep)

    call ref_blas_zgbmv_x(blas_no_trans, m, n, kl, ku, alpha, a, lda, &
                          x, 1, beta, y_ref, 1)
    call target_blas_zgbmv_x(blas_no_trans, m, n, kl, ku, alpha, a, lda, &
                             x, 1, beta, y_got, 1)

    err = max_rel_err_vec_z(y_got, y_ref)
    tol = 16.0_ep * (4.0_ep * real(m, ep) * real(kl + ku + 1, ep)) * target_eps
    write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm=', m, ' n=', n, ' kl=', kl, ' ku=', ku
    call report_case(trim(label), err, tol)

    call report_finalize()
end program test_blas_zgbmv_x
