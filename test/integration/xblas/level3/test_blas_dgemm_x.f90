program test_blas_dgemm_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: max_rel_err_mat
    use test_data,      only: gen_matrix_quad
    use target_xblas,   only: target_name, target_eps, target_blas_dgemm_x, &
                              blas_no_trans
    use ref_quad_xblas, only: ref_blas_dgemm_x
    implicit none

    integer, parameter :: cases(3, 3) = reshape([10,10,10, 32,24,16, 96,80,64], [3,3])
    integer :: c, m, n, k
    real(ep), allocatable :: a(:,:), b(:,:), c_ref(:,:), c_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call report_init('blas_dgemm_x', target_name)

    do c = 1, size(cases, 2)
        m = cases(1, c); n = cases(2, c); k = cases(3, c)
        call gen_matrix_quad(m, k, a, seed = 100 + 13 * c)
        call gen_matrix_quad(k, n, b, seed = 200 + 13 * c)
        call gen_matrix_quad(m, n, c_ref, seed = 300 + 13 * c)
        allocate(c_got(m, n))
        c_got = c_ref
        alpha = 1.7_ep
        beta  = -0.3_ep

        call ref_blas_dgemm_x(blas_no_trans, blas_no_trans, m, n, k, &
                              alpha, a, m, b, k, beta, c_ref, m)
        call target_blas_dgemm_x(blas_no_trans, blas_no_trans, m, n, k, &
                                 alpha, a, m, b, k, beta, c_got, m)

        err = max_rel_err_mat(c_got, c_ref)
        tol = 16.0_ep * (2.0_ep * real(m, ep) * real(n, ep) * real(k, ep)) * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ' n=', n, ' k=', k
        call report_case(trim(label), err, tol)
        deallocate(a, b, c_ref, c_got)
    end do

    call report_finalize()
end program test_blas_dgemm_x
