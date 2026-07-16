program test_blas_dgemv_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_vec
    use test_data,      only: gen_matrix_quad, gen_vector_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dgemv_x, &
                              blas_no_trans
    use ref_quad_xblas, only: ref_blas_dgemv_x
    implicit none

    integer, parameter :: cases(2, 3) = reshape([10,10, 64,48, 200,150], [2,3])
    integer :: c, m, n
    real(ep), allocatable :: a(:,:), x(:), y_ref(:), y_got(:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('blas_dgemv_x', target_name)

    do c = 1, size(cases, 2)
        m = cases(1, c); n = cases(2, c)
        call gen_matrix_quad(m, n, a, seed = 100 + 13 * c)
        call gen_vector_quad(n, x,    seed = 200 + 13 * c)
        call gen_vector_quad(m, y_ref, seed = 300 + 13 * c)
        allocate(y_got(m))
        y_got = y_ref
        alpha = 1.7_ep
        beta  = -0.3_ep

        call ref_blas_dgemv_x(blas_no_trans, m, n, alpha, a, m, x, 1, &
                              beta, y_ref, 1)
        call target_blas_dgemv_x(blas_no_trans, m, n, alpha, a, m, x, 1, &
                                 beta, y_got, 1)

        err = max_rel_err_vec(y_got, y_ref)
        tol = 16.0_ep * (2.0_ep * real(m, ep) * real(n, ep)) * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ' n=', n
        call report_case(trim(label), err, tol)
        deallocate(a, x, y_ref, y_got)
    end do

    call report_finalize()
end program test_blas_dgemv_x
