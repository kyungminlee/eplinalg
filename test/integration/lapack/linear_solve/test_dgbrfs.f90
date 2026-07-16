program test_dgbrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dgbtrf, target_dgbtrs, target_dgbrfs
    use ref_quad_lapack, only: dgbtrf, dgbtrs, dgbrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kl = 2, ku = 3, nrhs = 3
    integer :: i, n, ldab, ldafb, info, j, k
    integer, allocatable :: ipiv_ref(:), ipiv_got(:), iwork(:)
    real(ep), allocatable :: Adense(:,:), AB0(:,:), AFB_ref(:,:), AFB_got(:,:)
    real(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dgbrfs', target_name)
    do i = 1, size(ns)
        n = ns(i); ldab = kl + ku + 1; ldafb = 2*kl + ku + 1
        call gen_matrix_quad(n, n, Adense, seed = 138001 + 47 * i)
        do j = 1, n; Adense(j, j) = Adense(j, j) + real(2*n, ep); end do
        allocate(AB0(ldab, n)); AB0 = 0.0_ep
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AB0(ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        allocate(AFB_ref(ldafb, n), AFB_got(ldafb, n)); AFB_ref = 0.0_ep
        do j = 1, n
            do k = max(1, j-ku), min(n, j+kl)
                AFB_ref(kl + ku + 1 + k - j, j) = Adense(k, j)
            end do
        end do
        AFB_got = AFB_ref
        call gen_matrix_quad(n, nrhs, B0, seed = 138011 + 47 * i)
        allocate(ipiv_ref(n), ipiv_got(n), X_ref(n, nrhs), X_got(n, nrhs))
        allocate(work(3*n), iwork(n), ferr(nrhs), berr(nrhs))
        call dgbtrf(n, n, kl, ku, AFB_ref, ldafb, ipiv_ref, info)
        call target_dgbtrf(n, n, kl, ku, AFB_got, ldafb, ipiv_got, info)
        X_ref = B0; X_got = B0
        call dgbtrs('N', n, kl, ku, nrhs, AFB_ref, ldafb, ipiv_ref, X_ref, n, info)
        call target_dgbtrs('N', n, kl, ku, nrhs, AFB_got, ldafb, ipiv_got, X_got, n, info)
        call dgbrfs('N', n, kl, ku, nrhs, AB0, ldab, AFB_ref, ldafb, ipiv_ref, &
                    B0, n, X_ref, n, ferr, berr, work, iwork, info)
        call target_dgbrfs('N', n, kl, ku, nrhs, AB0, ldab, AFB_got, ldafb, ipiv_got, &
                           B0, n, X_got, n, ferr, berr, info)
        err = max_rel_err_mat(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Adense, AB0, AFB_ref, AFB_got, B0, ipiv_ref, ipiv_got, &
                   X_ref, X_got, work, iwork, ferr, berr)
    end do
    call report_finalize()
end program test_dgbrfs
