program test_dtrsm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat
    use test_data,     only: gen_matrix_quad
    use target_blas,   only: target_name, target_eps, target_dtrsm
    use ref_quad_blas, only: dtrsm
    implicit none

    ! Sweep SIDE × UPLO × TRANS × DIAG = 16 combos. Right-side and
    ! lower-triangle solves share zero code with the left-upper variant.
    integer, parameter :: ms(*) = [4, 32, 100, 64, 4, 32, 100, 64, 4, 32, 100, 64, 4, 32, 100, 64]
    integer, parameter :: ns(*) = [5, 40,  80, 48, 5, 40,  80, 48, 5, 40,  80, 48, 5, 40,  80, 48]
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R']
    character(len=1), parameter :: uplos(*)   = ['U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T', 'N', 'N', 'T', 'T', 'N', 'N', 'T', 'T', 'N', 'N', 'T', 'T']
    character(len=1), parameter :: diags(*)   = ['N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U']
    integer :: i, m, n, j, na, lda, ldb
    real(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=48) :: label

    call report_init('dtrsm', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        ! For SIDE='R' the triangular factor is n×n, not m×m.
        na = merge(m, n, sides(i) == 'L')
        call gen_matrix_quad(na, na, A, seed = 931 + 23 * i)
        ! Beef up the diagonal so the triangle is well-conditioned.
        do j = 1, na
            A(j, j) = A(j, j) + real(na + 2, ep)
        end do
        call gen_matrix_quad(m, n, B0, seed = 941 + 23 * i)
        alpha = real(0.7_ep, ep)
        allocate(B_ref(m, n), B_got(m, n))
        B_ref = B0
        B_got = B0
        lda = na
        ldb = m
        call dtrsm(sides(i), uplos(i), transes(i), diags(i), m, n, alpha, &
                   A, lda, B_ref, ldb)
        call target_dtrsm(sides(i), uplos(i), transes(i), diags(i), m, n, &
                          alpha, A, lda, B_got, ldb)
        err = max_rel_err_mat(B_got, B_ref)
        tol = 16.0_ep * 4.0_ep * real(na, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,a,a,i0,a,i0)') &
            'side=', sides(i), ',uplo=', uplos(i), ',trans=', transes(i), &
            ',diag=', diags(i), ',m=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B_ref, B_got)
    end do
    call report_finalize()
end program test_dtrsm
