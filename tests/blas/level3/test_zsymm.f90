program test_zsymm
    use prec_kinds,    only: ep
    use prec_report,   only: report_init, report_case, report_finalize
    use compare,       only: max_rel_err_mat_z
    use test_data,     only: gen_matrix_complex
    use target_blas,   only: target_name, target_eps, target_zsymm
    use ref_quad_blas, only: zsymm
    implicit none

    integer, parameter :: cases(*)              = [16, 64, 32, 48]
    ! Sweep SIDE × UPLO = 4 combos.
    character(len=1), parameter :: sides(*)    = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: uplos(*)    = ['U', 'L', 'U', 'L']
    integer :: i, m, n, an
    complex(ep), allocatable :: A(:,:), B(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=64) :: label

    call report_init('zsymm', target_name)
    do i = 1, size(cases)
        n = cases(i); m = cases(i)
        if (sides(i) == 'L') then; an = m; else; an = n; end if
        call gen_matrix_complex(an, an, A, seed = 2901 + 17 * i)
        call gen_matrix_complex(m,  n,  B, seed = 2911 + 17 * i)
        call gen_matrix_complex(m,  n,  C0, seed = 2921 + 17 * i)
        alpha = cmplx(0.6_ep, 0.2_ep, ep)
        beta  = cmplx(0.4_ep, -0.1_ep, ep)
        allocate(C_ref(m,n), C_got(m,n))
        C_ref = C0; C_got = C0
        call zsymm(sides(i), uplos(i), m, n, alpha, A, an, B, m, beta, C_ref, m)
        call target_zsymm(sides(i), uplos(i), m, n, alpha, A, an, B, m, beta, C_got, m)
        err = max_rel_err_mat_z(C_got, C_ref)
        tol = 64.0_ep * real(an, ep) * target_eps
        write(label, '(a,a,a,a,a,i0)') 'side=', sides(i), &
            ',uplo=', uplos(i), ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(C_ref, C_got)
    end do
    call report_finalize()
end program test_zsymm
