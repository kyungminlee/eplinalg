program test_dtbrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtbtrs, target_dtbrfs
    use ref_quad_lapack, only: dtbtrs, dtbrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kd = 3, nrhs = 3
    integer :: i, n, ldab, info, j, k
    real(ep), allocatable :: A(:,:), AB(:,:), B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: ferr_ref(:), berr_ref(:), ferr_got(:), berr_got(:)
    integer, allocatable :: iwork(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtbrfs', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kd + 1
        call gen_matrix_quad(n, n, A, seed = 156001 + 47 * i)
        do j = 1, n; A(j, j) = A(j, j) + real(2*n, ep); end do
        allocate(AB(ldab, n)); AB = 0.0_ep
        do j = 1, n
            do k = max(1, j-kd), j
                AB(kd + 1 + k - j, j) = A(k, j)
            end do
        end do
        call gen_matrix_quad(n, nrhs, B0, seed = 156011 + 47 * i)
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(3*n), iwork(n))
        allocate(ferr_ref(nrhs), berr_ref(nrhs), ferr_got(nrhs), berr_got(nrhs))
        X_ref = B0; X_got = B0
        call dtbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, X_ref, n, info)
        call target_dtbtrs('U', 'N', 'N', n, kd, nrhs, AB, ldab, X_got, n, info)
        call dtbrfs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B0, n, X_ref, n, &
                    ferr_ref, berr_ref, work, iwork, info)
        call target_dtbrfs('U', 'N', 'N', n, kd, nrhs, AB, ldab, B0, n, X_ref, n, &
                           ferr_got, berr_got, info)
        err = maxval(abs(berr_got - berr_ref))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, AB, B0, X_ref, X_got, work, iwork, ferr_ref, berr_ref, ferr_got, berr_got)
    end do
    call report_finalize()
end program test_dtbrfs
