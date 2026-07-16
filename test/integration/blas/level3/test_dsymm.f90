program test_dsymm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dsymm
    use ref_quad_blas, only: dsymm
    implicit none

    integer, parameter :: ms(*) = [4, 32, 100, 64]
    integer, parameter :: ns(*) = [5, 40,  80, 48]
    ! Cycle (SIDE, UPLO) so the four shapes touch every code path:
    ! left-upper (default), right-upper, left-lower, right-lower.
    character(len=1), parameter :: sides(*) = ['L', 'R', 'L', 'R']
    character(len=1), parameter :: uplos(*) = ['U', 'U', 'L', 'L']
    integer :: i, m, n, na, lda
    real(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call report_init('dsymm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! For SIDE='R' the symmetric factor is n×n, not m×m.
        na = merge(m, n, sides(i) == 'L')
        call gen_matrix_quad(na, na, A, seed = 831 + 23 * i)
        call gen_matrix_quad(m, n, B,  seed = 841 + 23 * i)
        call gen_matrix_quad(m, n, C0, seed = 851 + 23 * i)
        alpha = real(0.6_ep, ep)
        beta  = real(0.4_ep, ep)
        allocate(C_ref(m, n), C_got(m, n))
        C_ref = C0
        C_got = C0
        lda = na
        call dsymm(sides(i), uplos(i), m, n, alpha, A, lda, B, m, &
                   beta, C_ref, m)
        call target_dsymm(sides(i), uplos(i), m, n, alpha, A, lda, B, m, &
                          beta, C_got, m)
        err = max_rel_err_mat(C_got, C_ref)
        tol = 16.0_ep * 2.0_ep * real(na, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') &
            'side=', sides(i), ',uplo=', uplos(i), ',m=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, C0, C_ref, C_got)
    end do
    call report_finalize()
end program test_dsymm
