program test_zhemm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zhemm
    use ref_quad_blas, only: zhemm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 80, 64]
    integer, parameter :: ns(*) = [5, 40, 60, 48]
    ! Cycle (SIDE, UPLO) so the four shapes touch every Hermitian path.
    ! On SIDE='R' the Hermitian factor is n×n; the diagonal-conjugate
    ! handling differs between the upper and lower halves.
    character(len=1), parameter :: sides(*) = ['L', 'R', 'L', 'R']
    character(len=1), parameter :: uplos(*) = ['U', 'U', 'L', 'L']
    integer :: i, m, n, na, lda
    complex(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhemm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        na = merge(m, n, sides(i) == 'L')
        call gen_matrix_complex(na, na, A, seed = 1031 + 29 * i)
        call gen_matrix_complex(m, n, B,  seed = 1041 + 29 * i)
        call gen_matrix_complex(m, n, C0, seed = 1051 + 29 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        lda = na
        call zhemm(sides(i), uplos(i), m, n, alpha, A, lda, B, m, &
                   beta, C_ref, m)
        call target_zhemm(sides(i), uplos(i), m, n, alpha, A, lda, B, m, &
                          beta, C_got, m)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 32.0_ep * 8.0_ep * real(na, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') &
            'side=', sides(i), ',uplo=', uplos(i), ',m=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, C0, C_ref, C_got)
    end do
    call report_finalize()
end program test_zhemm
