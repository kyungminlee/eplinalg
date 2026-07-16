program test_zgtrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_vector_complex, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zgttrf, target_zgttrs, target_zgtrfs
    use ref_quad_lapack, only: zgttrf, zgttrs, zgtrfs
    implicit none

    integer, parameter :: ns(*) = [16, 64, 128]
    integer, parameter :: nrhs = 2
    integer :: i, n, info, j
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: dl0(:), d0(:), du0(:), B0(:,:)
    complex(ep), allocatable :: dlf_ref(:), df_ref(:), duf_ref(:), du2_ref(:)
    complex(ep), allocatable :: dlf_got(:), df_got(:), duf_got(:), du2_got(:)
    complex(ep), allocatable :: X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgtrfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_vector_complex(n-1, dl0, seed = 143001 + 47 * i)
        call gen_vector_complex(n,   d0,  seed = 143011 + 47 * i)
        call gen_vector_complex(n-1, du0, seed = 143021 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 143031 + 47 * i)
        do j = 1, n; d0(j) = d0(j) + cmplx(real(4, ep), 0.0_ep, ep); end do
        allocate(dlf_ref(n-1), df_ref(n), duf_ref(n-1), du2_ref(n-2))
        allocate(dlf_got(n-1), df_got(n), duf_got(n-1), du2_got(n-2))
        allocate(ipiv_ref(n), ipiv_got(n), X_ref(n, nrhs), X_got(n, nrhs))
        allocate(work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        dlf_ref = dl0; df_ref = d0; duf_ref = du0
        dlf_got = dl0; df_got = d0; duf_got = du0
        call zgttrf(n, dlf_ref, df_ref, duf_ref, du2_ref, ipiv_ref, info)
        call target_zgttrf(n, dlf_got, df_got, duf_got, du2_got, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zgttrs('N', n, nrhs, dlf_ref, df_ref, duf_ref, du2_ref, ipiv_ref, X_ref, n, info)
        call target_zgttrs('N', n, nrhs, dlf_got, df_got, duf_got, du2_got, ipiv_got, X_got, n, info)
        call zgtrfs('N', n, nrhs, dl0, d0, du0, dlf_ref, df_ref, duf_ref, du2_ref, ipiv_ref, &
                    B0, n, X_ref, n, ferr, berr, work, rwork, info)
        call target_zgtrfs('N', n, nrhs, dl0, d0, du0, dlf_got, df_got, duf_got, du2_got, ipiv_got, &
                           B0, n, X_got, n, ferr, berr, info)
        err = max_rel_err_mat_z(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(dl0, d0, du0, B0, dlf_ref, df_ref, duf_ref, du2_ref)
        deallocate(dlf_got, df_got, duf_got, du2_got, ipiv_ref, ipiv_got)
        deallocate(X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zgtrfs
