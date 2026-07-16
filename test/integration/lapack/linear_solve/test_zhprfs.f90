program test_zhprfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hermitian_matrix_quad, gen_matrix_complex, &
                                pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhptrf, target_zhptrs, target_zhprfs
    use ref_quad_lapack, only: zhptrf, zhptrs, zhprfs
    implicit none
    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, np, j
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), AP0(:), AP_ref(:), AP_got(:), B0(:,:), X_ref(:,:), X_got(:,:)
    complex(ep), allocatable :: work(:)
    real(ep),    allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('zhprfs', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A0, seed = 75001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 75011 + 47 * i)
        allocate(AP0(np), AP_ref(np), AP_got(np), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        call pack_herm_packed_quad('U', n, A0, AP0)
        AP_ref = AP0; AP_got = AP0
        call zhptrf('U', n, AP_ref, ipiv_ref, info)
        call target_zhptrf('U', n, AP_got, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zhptrs('U', n, nrhs, AP_ref, ipiv_ref, X_ref, n, info)
        call target_zhptrs('U', n, nrhs, AP_got, ipiv_got, X_got, n, info)
        call zhprfs('U', n, nrhs, AP0, AP_ref, ipiv_ref, B0, n, X_ref, n, &
                    ferr, berr, work, rwork, info)
        call target_zhprfs('U', n, nrhs, AP0, AP_got, ipiv_got, B0, n, X_got, n, &
                           ferr, berr, info)
        err = max_rel_err_mat_z(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP0, AP_ref, AP_got, ipiv_ref, ipiv_got, X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zhprfs
