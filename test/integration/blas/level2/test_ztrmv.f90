program test_ztrmv
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_vec_z
    use test_data,     only: gen_matrix_complex, gen_vector_complex
    use target_blas,   only: target_name, target_eps, target_ztrmv
    use ref_quad_blas, only: ztrmv
    implicit none

    integer, parameter :: cases(*)             = [10, 50, 200, 10, 50, 200, 10, 50, 200, 10, 50, 200]
    ! Sweep UPLO × TRANS × DIAG = 12 combos. The 'C' branch (conjugate-
    ! transpose) is its own code path in ztrmv.
    character(len=1), parameter :: uplos(*)    = ['U', 'U', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transes(*)  = ['N', 'N', 'T', 'T', 'C', 'C', 'N', 'N', 'T', 'T', 'C', 'C']
    character(len=1), parameter :: diags(*)    = ['N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U']
    integer :: i, n
    complex(ep), allocatable :: A(:,:), x0(:), x_ref(:), x_got(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztrmv', target_name)
    do i = 1, size(cases)
        n = cases(i)
        call gen_matrix_complex(n, n, A,  seed = 2501 + 17 * i)
        call gen_vector_complex(n,     x0, seed = 2511 + 17 * i)
        allocate(x_ref(n), x_got(n))
        x_ref = x0; x_got = x0
        call ztrmv(uplos(i), transes(i), diags(i), n, A, n, x_ref, 1)
        call target_ztrmv(uplos(i), transes(i), diags(i), n, A, n, x_got, 1)
        err = max_rel_err_vec_z(x_got, x_ref)
        tol = 32.0_ep * real(n, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,i0)') 'uplo=', uplos(i), &
            ',trans=', transes(i), ',diag=', diags(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(x_ref, x_got)
    end do
    call report_finalize()
end program test_ztrmv
