! zuncsd: CS decomposition of a 2x2-partitioned unitary matrix (complex).
program test_zuncsd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zuncsd
    use ref_quad_lapack, only: zgeqrf, zungqr, zuncsd
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, p, q, info, lwork, j
    complex(ep), allocatable :: A(:,:), tau(:), work(:)
    complex(ep), allocatable :: X11r(:,:), X12r(:,:), X21r(:,:), X22r(:,:)
    complex(ep), allocatable :: X11g(:,:), X12g(:,:), X21g(:,:), X22g(:,:)
    complex(ep), allocatable :: U1(:,:), U2(:,:), V1t(:,:), V2t(:,:)
    real(ep),    allocatable :: theta_r(:), theta_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    integer, allocatable :: iwork(:)
    character(len=48) :: label

    call report_init('zuncsd', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = m/2
        call gen_matrix_complex(m, m, A, seed = 24151 + 79 * i)
        allocate(tau(m))
        call zgeqrf(m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(m, m, A, m, tau, work, lwork, info)
        deallocate(work)
        call zungqr(m, m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zungqr(m, m, m, A, m, tau, work, lwork, info)
        deallocate(work, tau)
        allocate(X11r(p, q), X12r(p, m-q), X21r(m-p, q), X22r(m-p, m-q))
        allocate(X11g(p, q), X12g(p, m-q), X21g(m-p, q), X22g(m-p, m-q))
        do j = 1, q;   X11r(:, j) = A(1:p, j);   X21r(:, j) = A(p+1:m, j);   end do
        do j = 1, m-q; X12r(:, j) = A(1:p, q+j); X22r(:, j) = A(p+1:m, q+j); end do
        X11g = X11r; X12g = X12r; X21g = X21r; X22g = X22r
        deallocate(A)
        allocate(U1(p, p), U2(m-p, m-p), V1t(q, q), V2t(m-q, m-q), &
                 theta_r(q), theta_g(q), iwork(m))
        block
            real(ep), allocatable :: rwork(:)
            real(ep) :: ropt(1)
            integer  :: lrwork
            call zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'D', m, p, q, &
                        X11r, p, X12r, p, X21r, m-p, X22r, m-p, theta_r, &
                        U1, p, U2, m-p, V1t, q, V2t, m-q, &
                        wopt, -1, ropt, -1, iwork, info)
            lwork  = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
            lrwork = max(1, int(ropt(1)));            allocate(rwork(lrwork))
            call zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'D', m, p, q, &
                        X11r, p, X12r, p, X21r, m-p, X22r, m-p, theta_r, &
                        U1, p, U2, m-p, V1t, q, V2t, m-q, &
                        work, lwork, rwork, lrwork, iwork, info)
            deallocate(work, rwork)
        end block
        call target_zuncsd('Y', 'Y', 'Y', 'Y', 'N', 'D', m, p, q, &
                           X11g, p, X12g, p, X21g, m-p, X22g, m-p, theta_g, &
                           U1, p, U2, m-p, V1t, q, V2t, m-q, info)
        err = max_rel_err_vec(theta_g, theta_r)
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X12r, X21r, X22r, X11g, X12g, X21g, X22g, &
                   U1, U2, V1t, V2t, theta_r, theta_g, iwork)
    end do
    call report_finalize()
end program test_zuncsd
