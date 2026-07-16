program test_zspr
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_complex_symmetric_quad, gen_vector_complex
    use target_lapack, only: target_name, target_eps, target_zspr
    use ref_quad_lapack, only: zspr
    implicit none
    integer, parameter :: ns(*) = [10, 50, 100]
    character(len=1), parameter :: uplos(2) = ['U', 'L']
    integer :: i, n, ju, np, j, k, idx
    complex(ep), allocatable :: A(:,:), AP0(:), AP_ref(:), AP_got(:), x(:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zspr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A, seed = 70301 + 19 * i)
        call gen_vector_complex(n, x, seed = 70311 + 19 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        np = n*(n+1)/2
        allocate(AP0(np), AP_ref(np), AP_got(np))
        do ju = 1, size(uplos)
            idx = 0
            if (uplos(ju) == 'U') then
                do j = 1, n; do k = 1, j; idx = idx + 1; AP0(idx) = A(k, j); end do; end do
            else
                do j = 1, n; do k = j, n; idx = idx + 1; AP0(idx) = A(k, j); end do; end do
            end if
            AP_ref = AP0; AP_got = AP0
            call zspr(uplos(ju), n, alpha, x, 1, AP_ref)
            call target_zspr(uplos(ju), n, alpha, x, 1, AP_got)
            err = max_rel_err_vec_z(AP_got, AP_ref)
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(ju), ',n=', n
            call report_case(trim(label), err, tol)
        end do
        deallocate(AP0, AP_ref, AP_got)
    end do
    call report_finalize()
end program test_zspr
