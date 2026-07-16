! dbbcsd: bidiagonal CS decomposition.
! Driven from a small synthetic theta/phi sequence with identity U/V.
program test_dbbcsd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use target_lapack,   only: target_name, target_eps, target_dbbcsd
    use ref_quad_lapack, only: dbbcsd
    implicit none

    integer, parameter :: ms(*) = [12, 24]
    integer :: i, m, p, q, info, j
    real(ep), allocatable :: theta_r(:), theta_g(:), phi_r(:), phi_g(:)
    real(ep), allocatable :: U1(:,:), U2(:,:), V1t(:,:), V2t(:,:)
    real(ep), allocatable :: B11d(:), B11e(:), B12d(:), B12e(:)
    real(ep), allocatable :: B21d(:), B21e(:), B22d(:), B22e(:), work(:)
    real(ep) :: wopt(1), err, tol
    integer :: lwork
    character(len=48) :: label

    call report_init('dbbcsd', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = m/2
        allocate(theta_r(q), theta_g(q), phi_r(q-1), phi_g(q-1))
        allocate(U1(p, p), U2(m-p, m-p), V1t(q, q), V2t(m-q, m-q))
        allocate(B11d(q), B11e(q), B12d(q), B12e(q), B21d(q), B21e(q), B22d(q), B22e(q))
        do j = 1, q
            theta_r(j) = 0.1_ep + 0.7_ep * real(j-1, ep) / real(q, ep)
        end do
        theta_g = theta_r
        phi_r = 0.05_ep; phi_g = 0.05_ep
        U1 = 0.0_ep; U2 = 0.0_ep; V1t = 0.0_ep; V2t = 0.0_ep
        do j = 1, p; U1(j,j) = 1.0_ep; end do
        do j = 1, m-p; U2(j,j) = 1.0_ep; end do
        do j = 1, q; V1t(j,j) = 1.0_ep; end do
        do j = 1, m-q; V2t(j,j) = 1.0_ep; end do
        call dbbcsd('Y', 'Y', 'Y', 'Y', 'N', m, p, q, theta_r, phi_r, &
                    U1, p, U2, m-p, V1t, q, V2t, m-q, &
                    B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dbbcsd('Y', 'Y', 'Y', 'Y', 'N', m, p, q, theta_r, phi_r, &
                    U1, p, U2, m-p, V1t, q, V2t, m-q, &
                    B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, work, lwork, info)
        deallocate(work)
        deallocate(U1, U2, V1t, V2t, B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e)
        ! Reset for target call
        allocate(U1(p, p), U2(m-p, m-p), V1t(q, q), V2t(m-q, m-q))
        allocate(B11d(q), B11e(q), B12d(q), B12e(q), B21d(q), B21e(q), B22d(q), B22e(q))
        U1 = 0.0_ep; U2 = 0.0_ep; V1t = 0.0_ep; V2t = 0.0_ep
        do j = 1, p; U1(j,j) = 1.0_ep; end do
        do j = 1, m-p; U2(j,j) = 1.0_ep; end do
        do j = 1, q; V1t(j,j) = 1.0_ep; end do
        do j = 1, m-q; V2t(j,j) = 1.0_ep; end do
        call target_dbbcsd('Y', 'Y', 'Y', 'Y', 'N', m, p, q, theta_g, phi_g, &
                           U1, p, U2, m-p, V1t, q, V2t, m-q, &
                           B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e, info)
        err = max_rel_err_vec(theta_g, theta_r)
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(theta_r, theta_g, phi_r, phi_g, U1, U2, V1t, V2t, &
                   B11d, B11e, B12d, B12e, B21d, B21e, B22d, B22e)
    end do
    call report_finalize()
end program test_dbbcsd
