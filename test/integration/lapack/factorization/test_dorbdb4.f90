! dorbdb4: variant of orbdb for the case M-Q ≤ min(P, M-P, Q).
program test_dorbdb4
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorbdb4
    use ref_quad_lapack, only: dgeqrf, dorgqr, dorbdb4
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, p, q, info, lwork
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: X11r(:,:), X21r(:,:), X11g(:,:), X21g(:,:)
    real(ep), allocatable :: theta_r(:), theta_g(:), phi_r(:), phi_g(:)
    real(ep), allocatable :: tp1(:), tp2(:), tq1(:), phantom(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorbdb4', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = 3*m/4  ! M-Q ≤ P, M-Q ≤ M-P, M-Q ≤ Q
        call gen_matrix_quad(m, m, A, seed = 24601 + 79 * i)
        allocate(tau(m))
        call dgeqrf(m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(m, m, A, m, tau, work, lwork, info)
        deallocate(work)
        call dorgqr(m, m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorgqr(m, m, m, A, m, tau, work, lwork, info)
        deallocate(work, tau)
        allocate(X11r(p, q), X21r(m-p, q), X11g(p, q), X21g(m-p, q))
        X11r = A(1:p, 1:q); X21r = A(p+1:m, 1:q)
        X11g = X11r; X21g = X21r
        deallocate(A)
        allocate(theta_r(q), theta_g(q), phi_r(max(1,q-1)), phi_g(max(1,q-1)), &
                 tp1(max(1,m-q)), tp2(max(1,m-q)), tq1(q), phantom(m))
        call dorbdb4(m, p, q, X11r, p, X21r, m-p, theta_r, phi_r, &
                     tp1, tp2, tq1, phantom, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorbdb4(m, p, q, X11r, p, X21r, m-p, theta_r, phi_r, &
                     tp1, tp2, tq1, phantom, work, lwork, info)
        deallocate(work)
        call target_dorbdb4(m, p, q, X11g, p, X21g, m-p, theta_g, phi_g, &
                            tp1, tp2, tq1, phantom, info)
        ! Only theta(1:m-q) is filled by orbdb4 (loop runs 1..M-Q).
        err = max_rel_err_vec(theta_g(1:m-q), theta_r(1:m-q))
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X21r, X11g, X21g, theta_r, theta_g, phi_r, phi_g, tp1, tp2, tq1, phantom)
    end do
    call report_finalize()
end program test_dorbdb4
