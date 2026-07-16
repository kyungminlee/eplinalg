program test_dtrsv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dtrsv
    use ref_quad_blas, only: dtrsv
    implicit none

    integer, parameter :: cases(*) = [10, 50, 200, 64, 10, 50, 200, 64]
    ! Sweep the full UPLO × TRANS × DIAG product (8 combos). The diagonal
    ! is beefed up below so back-/forward-substitution stays well-
    ! conditioned for every shape; DIAG='U' ignores the stored diagonal.
    character(len=1), parameter :: uplos(*)   = ['U', 'U', 'U', 'U', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T', 'N', 'N', 'T', 'T']
    character(len=1), parameter :: diags(*)   = ['N', 'U', 'N', 'U', 'N', 'U', 'N', 'U']
    integer :: i, n, j
    real(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtrsv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_quad(n, n, A, seed = 641 + 17 * i)
        ! Beef up the diagonal so the triangle is well-conditioned for
        ! back- / forward-substitution. (For DIAG='U' the routine
        ! ignores diagonal entries, but we set them anyway for safety.)
        do j = 1, n
            A(j, j) = A(j, j) + real(n + 2, ep)
        end do
        call gen_vector_quad(n, x0, seed = 651 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0
        x_got = x0
        call dtrsv(uplos(i), transes(i), diags(i), n, A, n, x_ref, 1)
        call target_dtrsv(uplos(i), transes(i), diags(i), n, A, n, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 16.0_ep * 4.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',diag=', diags(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_dtrsv
