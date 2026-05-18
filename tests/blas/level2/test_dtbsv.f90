program test_dtbsv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec
    use test_data,     only: gen_matrix_quad, gen_vector_quad
    use target_blas,   only: target_name, target_eps, target_dtbsv
    use ref_quad_blas, only: dtbsv
    implicit none

    ! DIAG='U' for every case so the (implicit-1) diagonal keeps the
    ! solve well-conditioned regardless of the random off-diagonals.
    integer, parameter :: cases(*)              = [20, 100, 64, 50]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'N', 'T', 'T']
    integer,          parameter :: ks(*)       = [3, 5, 0, 2]
    integer :: i, n, k, lda
    real(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtbsv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        k   = ks(i)
        lda = k + 1
        call gen_matrix_quad(lda, n, A,  seed = 1101 + 17 * i)
        ! Shrink off-diagonals so the unit-diag back-solve stays bounded.
        A = 0.1_ep * A
        call gen_vector_quad(n,      x0, seed = 1111 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call dtbsv(uplos(i), transes(i), 'U', n, k, A, lda, x_ref, 1)
        call target_dtbsv(uplos(i), transes(i), 'U', n, k, A, lda, x_got, 1)
        err = max_rel_err_vec(x_got, x_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, x_ref, x_got)
    end do
    call report_finalize()
end program test_dtbsv
