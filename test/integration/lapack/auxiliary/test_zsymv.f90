program test_zsymv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_complex_symmetric_quad, gen_vector_complex
    use target_lapack, only: target_name, target_eps, target_zsymv
    use ref_quad_lapack, only: zsymv
    implicit none
    integer, parameter :: ns(*) = [10, 50, 100]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, ju
    complex(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zsymv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A, seed = 70001 + 19 * i)
        call gen_vector_complex(n, x,  seed = 70011 + 19 * i)
        call gen_vector_complex(n, y0, seed = 70021 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep); beta = cmplx(0.4_ep, -0.1_ep, ep)
        do ju = 1, size(uplos)
            allocate(y_ref(n), y_got(n)); y_ref = y0; y_got = y0
            call zsymv(uplos(ju), n, alpha, A, n, x, 1, beta, y_ref, 1)
            call target_zsymv(uplos(ju), n, alpha, A, n, x, 1, beta, y_got, 1)
            err = max_rel_err_vec_z(y_got, y_ref)
            tol = 32.0_ep * 4.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(y_ref, y_got)
        end do
    end do
    call report_finalize()
end program test_zsymv
