! dorcsd2by1: 2-by-1 CS decomposition of [X11; X21] orthonormal block.
program test_dorcsd2by1
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorcsd2by1
    use ref_quad_lapack, only: dgeqrf, dorgqr, dorcsd2by1
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, p, q, info, lwork
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: X11r(:,:), X21r(:,:), X11g(:,:), X21g(:,:)
    real(ep), allocatable :: U1(:,:), U2(:,:), V1t(:,:)
    real(ep), allocatable :: theta_r(:), theta_g(:)
    real(ep) :: wopt(1), err, tol
    integer, allocatable :: iwork(:)
    character(len=48) :: label

    call report_init('dorcsd2by1', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = m/2
        call gen_matrix_quad(m, m, A, seed = 24201 + 79 * i)
        allocate(tau(m))
        call dgeqrf(m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(m, m, A, m, tau, work, lwork, info)
        deallocate(work)
        call dorgqr(m, m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorgqr(m, m, m, A, m, tau, work, lwork, info)
        deallocate(work, tau)
        ! Slice first q columns of orthogonal m×m Q into [X11; X21].
        allocate(X11r(p, q), X21r(m-p, q), X11g(p, q), X21g(m-p, q))
        X11r = A(1:p, 1:q); X21r = A(p+1:m, 1:q)
        X11g = X11r; X21g = X21r
        deallocate(A)
        allocate(U1(p, p), U2(m-p, m-p), V1t(q, q), &
                 theta_r(q), theta_g(q), iwork(m-q))
        call dorcsd2by1('Y', 'Y', 'Y', m, p, q, X11r, p, X21r, m-p, &
                        theta_r, U1, p, U2, m-p, V1t, q, &
                        wopt, -1, iwork, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorcsd2by1('Y', 'Y', 'Y', m, p, q, X11r, p, X21r, m-p, &
                        theta_r, U1, p, U2, m-p, V1t, q, &
                        work, lwork, iwork, info)
        deallocate(work)
        call target_dorcsd2by1('Y', 'Y', 'Y', m, p, q, X11g, p, X21g, m-p, &
                               theta_g, U1, p, U2, m-p, V1t, q, info)
        err = max_rel_err_vec(theta_g, theta_r)
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X21r, X11g, X21g, U1, U2, V1t, theta_r, theta_g, iwork)
    end do
    call report_finalize()
end program test_dorcsd2by1
