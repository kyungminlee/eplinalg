! dtrrfs computes FERR and BERR for a triangular solve. Compare BERR.
program test_dtrrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrtrs, target_dtrrfs
    use ref_quad_lapack, only: dtrtrs, dtrrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, j
    real(ep), allocatable :: A(:,:), B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:)
    integer, allocatable :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtrrfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 152001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        call gen_matrix_quad(n, nrhs, B0, seed = 152011 + 47 * i)
        allocate(X_ref(n, nrhs), X_got(n, nrhs))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        allocate(work(3*n), iwork(n))
        X_ref = B0; X_got = B0
        call dtrtrs('U', 'N', 'N', n, nrhs, A, n, X_ref, n, info)
        call target_dtrtrs('U', 'N', 'N', n, nrhs, A, n, X_got, n, info)
        call dtrrfs('U', 'N', 'N', n, nrhs, A, n, B0, n, X_ref, n, &
                    ferr_ref, berr_ref, work, iwork, info)
        call target_dtrrfs('U', 'N', 'N', n, nrhs, A, n, B0, n, X_ref, n, &
                           ferr_got, berr_got, info)
        err = maxval(abs(berr_got - berr_ref))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, X_ref, X_got, ferr_ref, berr_ref, ferr_got, berr_got, work, iwork)
    end do
    call report_finalize()
end program test_dtrrfs
