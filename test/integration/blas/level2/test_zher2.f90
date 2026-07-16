program test_zher2
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_zher2
    use ref_quad_blas, only: zher2
    implicit none

    integer, parameter :: cases(*)            = [10, 50, 200]
    character(len=1), parameter :: uplos(*)  = ['U', 'L', 'U']
    integer :: i, n
    complex(ep), allocatable :: A0(:,:), x(:), y(:), A_ref(:,:), A_got(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zher2', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_complex(n, n, A0, seed = 1701 + 17 * i)
        call gen_vector_complex(n,     x, seed = 1711 + 17 * i)
        call gen_vector_complex(n,     y, seed = 1721 + 17 * i)
        alpha = cmplx(0.65_ep, -0.2_ep, ep)
        allocate(A_ref(n,n), A_got(n,n))
        A_ref = A0; A_got = A0
        call zher2(uplos(i), n, alpha, x, 1, y, 1, A_ref, n)
        call target_zher2(uplos(i), n, alpha, x, 1, y, 1, A_got, n)
        err = max_rel_err_mat_z(A_got, A_ref)
        tol = 16.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got)
    end do
    call report_finalize()
end program test_zher2
