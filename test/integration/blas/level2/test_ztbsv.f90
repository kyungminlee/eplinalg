program test_ztbsv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_ztbsv
    use ref_quad_blas, only: ztbsv
    implicit none

    ! Sweep UPLO × TRANS = 6 combos. DIAG='U' is pinned (implicit-1
    ! diagonal keeps the solve bounded regardless of off-diagonals).
    integer, parameter :: cases(*)              = [20, 100, 64, 50, 20, 100]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'N', 'T', 'T', 'C', 'C']
    integer,          parameter :: ks(*)       = [3, 5, 0, 2, 3, 5]
    integer :: i, n, k, lda
    complex(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztbsv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        k   = ks(i)
        lda = k + 1
        call gen_matrix_complex(lda, n, A,  seed = 2201 + 17 * i)
        A = 0.1_ep * A
        call gen_vector_complex(n,      x0, seed = 2211 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call ztbsv(uplos(i), transes(i), 'U', n, k, A, lda, x_ref, 1)
        call target_ztbsv(uplos(i), transes(i), 'U', n, k, A, lda, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, x_ref, x_got)
    end do
    call report_finalize()
end program test_ztbsv
