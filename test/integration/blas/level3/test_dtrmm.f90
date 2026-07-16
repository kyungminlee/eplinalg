program test_dtrmm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dtrmm
    use ref_quad_blas, only: dtrmm
    implicit none

    ! Sweep SIDE × UPLO × TRANS × DIAG = 16 combos.
    integer, parameter :: ms(*) = [4, 32, 100, 64, 4, 32, 100, 64, 4, 32, 100, 64, 4, 32, 100, 64]
    integer, parameter :: ns(*) = [5, 40,  80, 48, 5, 40,  80, 48, 5, 40,  80, 48, 5, 40,  80, 48]
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R']
    character(len=1), parameter :: uplos(*)   = ['U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T', 'N', 'N', 'T', 'T', 'N', 'N', 'T', 'T', 'N', 'N', 'T', 'T']
    character(len=1), parameter :: diags(*)   = ['N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U']
    integer :: i, m, n, na, lda
    real(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=48) :: label

    call report_init('dtrmm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! Triangular factor is m×m for SIDE='L', n×n for SIDE='R'.
        na = merge(m, n, sides(i) == 'L')
        call gen_matrix_quad(na, na, A, seed = 911 + 23 * i)
        call gen_matrix_quad(m, n, B0, seed = 921 + 23 * i)
        alpha = real(0.7_ep, ep)
        allocate(B_ref(m, n), B_got(m, n))
        B_ref = B0
        B_got = B0
        lda = na
        call dtrmm(sides(i), uplos(i), transes(i), diags(i), m, n, alpha, &
                   A, lda, B_ref, m)
        call target_dtrmm(sides(i), uplos(i), transes(i), diags(i), m, n, &
                          alpha, A, lda, B_got, m)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * 2.0_ep * real(na, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,a,a,i0,a,i0)') &
            'side=', sides(i), ',uplo=', uplos(i), ',trans=', transes(i), &
            ',diag=', diags(i), ',m=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, B_ref, B_got)
    end do
    call report_finalize()
end program test_dtrmm
