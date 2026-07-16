program test_zspmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_complex_symmetric_quad, gen_vector_complex
    use target_lapack, only: target_name, target_eps, target_zspmv
    use ref_quad_lapack, only: zspmv
    implicit none
    integer, parameter :: ns(*) = [10, 50, 100]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, ju, np, j, k, idx
    complex(ep), allocatable :: A(:,:), AP(:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zspmv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A, seed = 70201 + 19 * i)
        call gen_vector_complex(n, x,  seed = 70211 + 19 * i)
        call gen_vector_complex(n, y0, seed = 70221 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep); beta = cmplx(0.4_ep, -0.1_ep, ep)
        np = n*(n+1)/2
        allocate(AP(np))
        do ju = 1, size(uplos)
            idx = 0
            if (uplos(ju) == 'U') then
                do j = 1, n; do k = 1, j; idx = idx + 1; AP(idx) = A(k, j); end do; end do
            else
                do j = 1, n; do k = j, n; idx = idx + 1; AP(idx) = A(k, j); end do; end do
            end if
            allocate(y_ref(n), y_got(n)); y_ref = y0; y_got = y0
            call zspmv(uplos(ju), n, alpha, AP, x, 1, beta, y_ref, 1)
            call target_zspmv(uplos(ju), n, alpha, AP, x, 1, beta, y_got, 1)
            err = max_rel_err_vec_z(y_got, y_ref)
            tol = 32.0_ep * 4.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(y_ref, y_got)
        end do
        deallocate(AP)
    end do
    call report_finalize()
end program test_zspmv
