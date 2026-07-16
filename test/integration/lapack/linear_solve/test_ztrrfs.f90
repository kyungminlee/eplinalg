program test_ztrrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrtrs, target_ztrrfs
    use ref_quad_lapack, only: ztrtrs, ztrrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: nrhs = 3
    integer :: i, n, info, j
    complex(ep), allocatable :: A(:,:), B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: rwork(:), ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztrrfs', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 153001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + cmplx(real(2*n, ep), 0.0_ep, ep); end do
        call gen_matrix_complex(n, nrhs, B0, seed = 153011 + 47 * i)
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), rwork(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        X_ref = B0; X_got = B0
        call ztrtrs('U', 'N', 'N', n, nrhs, A, n, X_ref, n, info)
        call target_ztrtrs('U', 'N', 'N', n, nrhs, A, n, X_got, n, info)
        call ztrrfs('U', 'N', 'N', n, nrhs, A, n, B0, n, X_ref, n, &
                    ferr_ref, berr_ref, work, rwork, info)
        call target_ztrrfs('U', 'N', 'N', n, nrhs, A, n, B0, n, X_ref, n, &
                           ferr_got, berr_got, info)
        err = maxval(abs(berr_got - berr_ref))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B0, X_ref, X_got, work, rwork, ferr_ref, berr_ref, ferr_got, berr_got)
    end do
    call report_finalize()
end program test_ztrrfs
