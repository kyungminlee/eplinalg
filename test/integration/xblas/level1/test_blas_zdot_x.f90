program test_blas_zdot_x
    use prec_kinds,     only: ep
    use prec_report,    only: report_init, report_case, report_finalize
    use compare,        only: rel_err_scalar_z
    use test_data,      only: gen_vector_complex
    use target_xblas,   only: target_name, target_eps, target_blas_zdot_x, &
                              blas_no_conj => blas_no_conj, &
                              blas_conj_k  => blas_conj_k
    use ref_quad_xblas, only: ref_blas_zdot_x
    implicit none

    integer, parameter :: cases(3) = [10, 100, 1000]
    integer :: i, n, conj
    complex(ep), allocatable :: x(:), y(:)
    complex(ep) :: alpha, beta, ref_r, got_r
    real(ep) :: err, tol
    character(len=32) :: label

    call report_init('blas_zdot_x', target_name)

    do i = 1, size(cases)
        n = cases(i)
        call gen_vector_complex(n, x, seed = 42 + 7 * i)
        call gen_vector_complex(n, y, seed = 99 + 7 * i)
        alpha = (1.7_ep, -0.4_ep)
        beta  = (-0.3_ep, 0.6_ep)

        do conj = 0, 1
            ref_r = (0.5_ep, 0.25_ep)
            got_r = (0.5_ep, 0.25_ep)
            call ref_blas_zdot_x(merge(blas_conj_k, blas_no_conj, conj == 1), &
                                 n, alpha, x, 1, beta, y, 1, ref_r)
            call target_blas_zdot_x(merge(blas_conj_k, blas_no_conj, conj == 1), &
                                    n, alpha, x, 1, beta, y, 1, got_r)
            err = rel_err_scalar_z(got_r, ref_r)
            tol = 16.0_ep * (4.0_ep * real(n, ep)) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ' conj=', conj
            call report_case(trim(label), err, tol)
        end do
    end do

    call report_finalize()
end program test_blas_zdot_x
