program test_ztrmm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_ztrmm
    use ref_quad_blas, only: ztrmm
    implicit none

    ! Sweep SIDE × UPLO × TRANS × DIAG = 24 combos. 'C' (conjugate-
    ! transpose) is its own complex-only code path.
    integer, parameter :: cases(*)             = [16, 64, 32, 24, 16, 64, 32, 24, 16, 64, 32, 24, 16, 64, 32, 24, 16, 64, 32, 24, 16, 64, 32, 24]
    character(len=1), parameter :: sides(*)    = ['L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'L', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R', 'R']
    character(len=1), parameter :: uplos(*)    = ['U', 'U', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'L', 'L', 'U', 'U', 'U', 'U', 'U', 'U', 'L', 'L', 'L', 'L', 'L', 'L']
    character(len=1), parameter :: transas(*)  = ['N', 'N', 'T', 'T', 'C', 'C', 'N', 'N', 'T', 'T', 'C', 'C', 'N', 'N', 'T', 'T', 'C', 'C', 'N', 'N', 'T', 'T', 'C', 'C']
    character(len=1), parameter :: diags(*)    = ['N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U', 'N', 'U']
    integer :: i, m, n, an
    complex(ep), allocatable :: A(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('ztrmm', target_name)
    do i = 1, size(cases)
        m = cases(i); n = cases(i)
        if (sides(i) == 'L') then; an = m; else; an = n; end if
        call gen_matrix_complex(an, an, A, seed = 3301 + 17 * i)
        call gen_matrix_complex(m,  n,  B0, seed = 3311 + 17 * i)
        alpha = cmplx(0.7_ep, 0.3_ep, ep)
        allocate(B_ref(m,n), B_got(m,n))
        B_ref = B0; B_got = B0
        call ztrmm(sides(i), uplos(i), transas(i), diags(i), m, n, alpha, A, an, B_ref, m)
        call target_ztrmm(sides(i), uplos(i), transas(i), diags(i), m, n, alpha, A, an, B_got, m)
        err = max_rel_err_mat_z(B_got, B_ref)
        tol = 64.0_ep * real(an, ep) * target_eps
        write(label, '(a,a,a,a,a,a,a,a,a,i0)') 'side=', sides(i), &
            ',uplo=', uplos(i), ',trans=', transas(i), ',diag=', diags(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(B_ref, B_got)
    end do
    call report_finalize()
end program test_ztrmm
