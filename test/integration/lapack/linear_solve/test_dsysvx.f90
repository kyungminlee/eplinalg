program test_dsysvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, rel_err_scalar
    use test_data,       only: gen_symmetric_matrix_quad, gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsysvx
    use ref_quad_lapack, only: dsysvx
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, lwork, j
    real(ep), allocatable :: A(:,:), B0(:,:)
    real(ep), allocatable :: AF_ref(:,:), X_ref(:,:), AF_got(:,:), X_got(:,:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), work(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep) :: wopt(1), rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dsysvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A, seed = 420041 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        call gen_matrix_quad(n, nrhs, B0, seed = 420051 + 47 * i)
        allocate(AF_ref(n, n), AF_got(n, n), X_ref(n, nrhs), X_got(n, nrhs))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(ipiv_ref(n), ipiv_got(n), iwork(n))
        call dsysvx('N', 'U', n, nrhs, A, n, AF_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, rcond_ref, ferr_ref, berr_ref, wopt, -1, iwork, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsysvx('N', 'U', n, nrhs, A, n, AF_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, rcond_ref, ferr_ref, berr_ref, work, lwork, iwork, info)
        call target_dsysvx('N', 'U', n, nrhs, A, n, AF_got, n, ipiv_got, B0, n, &
                           X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AF_ref, AF_got, X_ref, X_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, ipiv_ref, ipiv_got, iwork, work)
    end do
    call report_finalize()
end program test_dsysvx
