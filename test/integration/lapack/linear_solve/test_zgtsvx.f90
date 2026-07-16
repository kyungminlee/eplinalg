program test_zgtsvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_vector_complex, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgtsvx
    use ref_quad_lapack, only: zgtsvx
    implicit none
    integer, parameter :: ns(*) = [16, 32]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, j
    complex(ep), allocatable :: dl(:), d(:), du(:), B0(:,:)
    complex(ep), allocatable :: dlf_r(:), df_r(:), duf_r(:), du2_r(:), X_r(:,:)
    complex(ep), allocatable :: dlf_g(:), df_g(:), duf_g(:), du2_g(:), X_g(:,:), work(:)
    real(ep),    allocatable :: rwork(:), ferr(:), berr(:)
    integer,     allocatable :: ipiv_r(:), ipiv_g(:)
    real(ep) :: rcond_r, rcond_g, err, tol
    character(len=48) :: label
    call report_init('zgtsvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_complex(n,   d,  seed = 91001 + 19 * i)
        call gen_vector_complex(n-1, dl, seed = 91011 + 19 * i)
        call gen_vector_complex(n-1, du, seed = 91021 + 19 * i)
        do j = 1, n; d(j) = d(j) + cmplx(4.0_ep, 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 91031 + 19 * i)
        allocate(dlf_r(max(1,n-1)), df_r(n), duf_r(max(1,n-1)), du2_r(max(1,n-2)), X_r(n, nrhs))
        allocate(dlf_g(max(1,n-1)), df_g(n), duf_g(max(1,n-1)), du2_g(max(1,n-2)), X_g(n, nrhs))
        allocate(ipiv_r(n), ipiv_g(n), work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        call zgtsvx('N', 'N', n, nrhs, dl, d, du, dlf_r, df_r, duf_r, du2_r, ipiv_r, &
                    B0, n, X_r, n, rcond_r, ferr, berr, work, rwork, info)
        call target_zgtsvx('N', 'N', n, nrhs, dl, d, du, dlf_g, df_g, duf_g, du2_g, ipiv_g, &
                           B0, n, X_g, n, rcond_g, ferr, berr, info)
        err = max(max_rel_err_mat_z(X_g, X_r), rel_err_scalar(rcond_g, rcond_r))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dlf_r, df_r, duf_r, du2_r, X_r, dlf_g, df_g, duf_g, du2_g, X_g, ipiv_r, ipiv_g, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zgtsvx
