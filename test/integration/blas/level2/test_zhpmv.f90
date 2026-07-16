program test_zhpmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zhpmv
    use ref_quad_blas, only: zhpmv
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n, aps, j, idx
    complex(ep), allocatable :: ap(:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zhpmv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        aps = n * (n + 1) / 2
        call gen_vector_complex(aps, ap, seed = 1501 + 17 * i)
        call gen_vector_complex(n,    x,  seed = 1511 + 17 * i)
        call gen_vector_complex(n,    y0, seed = 1521 + 17 * i)
        ! Hermitian-packed contract: diagonal entries must be real.
        if (uplos(i) == 'U' .or. uplos(i) == 'u') then
            idx = 0
            do j = 1, n
                idx = idx + j   ! AP(j*(j+1)/2) holds A(j,j)
                ap(idx) = cmplx(real(ap(idx), ep), 0.0_ep, ep)
            end do
        else
            idx = 1
            do j = 1, n
                ap(idx) = cmplx(real(ap(idx), ep), 0.0_ep, ep)
                idx = idx + (n - j + 1)
            end do
        end if
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0; y_got = y0
        call zhpmv(uplos(i), n, alpha, ap, x, 1, beta, y_ref, 1)
        call target_zhpmv(uplos(i), n, alpha, ap, x, 1, beta, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zhpmv
