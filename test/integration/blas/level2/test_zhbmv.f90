program test_zhbmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zhbmv
    use ref_quad_blas, only: zhbmv
    implicit none

    integer, parameter :: cases(*)            = [20, 100, 64]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer,          parameter :: ks(*)     = [3, 5, 0]
    integer :: i, n, k, lda, j
    complex(ep), allocatable :: A(:,:), x(:), y0(:), y_ref(:), y_got(:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zhbmv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        k   = ks(i)
        lda = k + 1
        call gen_matrix_complex(lda, n, A,  seed = 1401 + 17 * i)
        call gen_vector_complex(n,      x,  seed = 1411 + 17 * i)
        call gen_vector_complex(n,      y0, seed = 1421 + 17 * i)
        ! Hermitian-banded contract: zero imaginary part on the diagonal row.
        if (uplos(i) == 'U' .or. uplos(i) == 'u') then
            do j = 1, n
                A(k + 1, j) = cmplx(real(A(k + 1, j), ep), 0.0_ep, ep)
            end do
        else
            do j = 1, n
                A(1, j) = cmplx(real(A(1, j), ep), 0.0_ep, ep)
            end do
        end if
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(y_ref(n), y_got(n))
        y_ref = y0; y_got = y0
        call zhbmv(uplos(i), n, k, alpha, A, lda, x, 1, beta, y_ref, 1)
        call target_zhbmv(uplos(i), n, k, alpha, A, lda, x, 1, beta, y_got, 1)
        err = max_rel_err_vec_z(y_got, y_ref)
        tol = 32.0_ep * real(k + 1, ep) * target_eps
        write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(y_ref, y_got)
    end do
    call report_finalize()
end program test_zhbmv
