program test_zsysvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, rel_err_scalar
    use test_data,       only: gen_complex_symmetric_quad, gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zsysvx
    use ref_quad_lapack, only: zsysvx
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer, parameter :: nrhs  = 2
    integer :: i, n, info, lwork, j
    complex(ep), allocatable :: A(:,:), B0(:,:)
    complex(ep), allocatable :: AF_ref(:,:), X_ref(:,:), AF_got(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:), rwork(:)
    integer,  allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep) :: wopt(1)
    real(ep) :: rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('zsysvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_complex_symmetric_quad(n, A, seed = 420061 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 420071 + 47 * i)
        allocate(AF_ref(n, n), AF_got(n, n), X_ref(n, nrhs), X_got(n, nrhs))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(ipiv_ref(n), ipiv_got(n), rwork(n))
        call zsysvx('N', 'U', n, nrhs, A, n, AF_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, rcond_ref, ferr_ref, berr_ref, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zsysvx('N', 'U', n, nrhs, A, n, AF_ref, n, ipiv_ref, B0, n, &
                    X_ref, n, rcond_ref, ferr_ref, berr_ref, work, lwork, rwork, info)
        call target_zsysvx('N', 'U', n, nrhs, A, n, AF_got, n, ipiv_got, B0, n, &
                           X_got, n, rcond_got, ferr_got, berr_got, info)
        err = max(max_rel_err_mat_z(X_got, X_ref), rel_err_scalar(rcond_got, rcond_ref))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, AF_ref, AF_got, X_ref, X_got, &
                   ferr_ref, berr_ref, ferr_got, berr_got, ipiv_ref, ipiv_got, rwork, work)
    end do
    call report_finalize()
end program test_zsysvx
