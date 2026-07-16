program test_zhemv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zhemv
    use ref_quad_blas, only: zhemv
    implicit none

    integer, parameter :: cases(*) = [10, 50, 100]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n
    complex(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhemv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_complex(n, n, A, seed = 741 + 19 * i)
        call gen_vector_complex(n, x,  seed = 751 + 19 * i)
        call gen_vector_complex(n, y0, seed = 761 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0
        y_got = y0
        call zhemv(uplos(i), n, alpha, A, n, x, 1, beta, y_ref, 1)
        call target_zhemv(uplos(i), n, alpha, A, n, x, 1, beta, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * 4.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zhemv
