program test_zgerfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, &
                                target_zgetrf, target_zgetrs, target_zgerfs
    use ref_quad_lapack, only: zgetrf, zgetrs, zgerfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:)
    complex(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgerfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A0, seed = 135001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 135011 + 47 * i)
        allocate(A_ref(n,n), A_got(n,n), ipiv_ref(n), ipiv_got(n))
        allocate(X_ref(n, nrhs), X_got(n, nrhs))
        allocate(work(2*n), rwork(n), ferr(nrhs), berr(nrhs))
        A_ref = A0; A_got = A0
        call zgetrf(n, n, A_ref, n, ipiv_ref, info)
        call target_zgetrf(n, n, A_got, n, ipiv_got, info)
        X_ref = B0; X_got = B0
        call zgetrs('N', n, nrhs, A_ref, n, ipiv_ref, X_ref, n, info)
        call target_zgetrs('N', n, nrhs, A_got, n, ipiv_got, X_got, n, info)
        call zgerfs('N', n, nrhs, A0, n, A_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, ferr, berr, work, rwork, info)
        call target_zgerfs('N', n, nrhs, A0, n, A_got, n, ipiv_got, B0, n, &
                           X_got, n, ferr, berr, info)
        err = max_rel_err_mat_z(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, A_got, ipiv_ref, ipiv_got, X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zgerfs
