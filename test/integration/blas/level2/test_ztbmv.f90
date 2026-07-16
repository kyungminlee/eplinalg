program test_ztbmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_ztbmv
    use ref_quad_blas, only: ztbmv
    implicit none

    integer, parameter :: cases(*)              = [20, 100, 64, 50]
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'C', 'T', 'N']
    character(len=1), parameter :: diags(*)    = ['N', 'N', 'U', 'N']
    integer,          parameter :: ks(*)       = [3, 5, 0, 2]
    integer :: i, n, k, lda
    complex(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztbmv', target_name)
    do i = 1, size(cases)
        n   = cases(i)
        k   = ks(i)
        lda = k + 1
        call gen_matrix_complex(lda, n, A,  seed = 2101 + 17 * i)
        call gen_vector_complex(n,      x0, seed = 2111 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call ztbmv(uplos(i), transes(i), diags(i), n, k, A, lda, x_ref, 1)
        call target_ztbmv(uplos(i), transes(i), diags(i), n, k, A, lda, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 32.0_ep * real(k + 1, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',diag=', diags(i), ',n=', n, ',k=', k
        call report_case(trim(label), err, tol)
        deallocate(A, x_ref, x_got)
    end do
    call report_finalize()
end program test_ztbmv
