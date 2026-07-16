! zunbdb3: complex unbdb variant. M-P ≤ MIN(P, Q, M-Q).
program test_zunbdb3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zunbdb3
    use ref_quad_lapack, only: zgeqrf, zungqr, zunbdb3
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, p, q, info, lwork
    complex(ep), allocatable :: A(:,:), tau(:), work(:)
    complex(ep), allocatable :: X11r(:,:), X21r(:,:), X11g(:,:), X21g(:,:)
    real(ep), allocatable :: theta_r(:), theta_g(:), phi_r(:), phi_g(:)
    complex(ep), allocatable :: tp1(:), tp2(:), tq1(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunbdb3', target_name)
    do i = 1, size(ms)
        m = ms(i); p = 3*m/4; q = m/2
        call gen_matrix_complex(m, m, A, seed = 24551 + 79 * i)
        allocate(tau(m))
        call zgeqrf(m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(m, m, A, m, tau, work, lwork, info)
        deallocate(work)
        call zungqr(m, m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zungqr(m, m, m, A, m, tau, work, lwork, info)
        deallocate(work, tau)
        allocate(X11r(p, q), X21r(m-p, q), X11g(p, q), X21g(m-p, q))
        X11r = A(1:p, 1:q); X21r = A(p+1:m, 1:q)
        X11g = X11r; X21g = X21r
        deallocate(A)
        allocate(theta_r(q), theta_g(q), phi_r(max(1,q-1)), phi_g(max(1,q-1)), &
                 tp1(p), tp2(m-p), tq1(q))
        call zunbdb3(m, p, q, X11r, p, X21r, m-p, theta_r, phi_r, &
                     tp1, tp2, tq1, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zunbdb3(m, p, q, X11r, p, X21r, m-p, theta_r, phi_r, &
                     tp1, tp2, tq1, work, lwork, info)
        deallocate(work)
        call target_zunbdb3(m, p, q, X11g, p, X21g, m-p, theta_g, phi_g, &
                            tp1, tp2, tq1, info)
        ! Only theta(1:m-p) is filled by unbdb3 (loop runs 1..M-P).
        err = max_rel_err_vec(theta_g(1:m-p), theta_r(1:m-p))
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X21r, X11g, X21g, theta_r, theta_g, phi_r, phi_g, tp1, tp2, tq1)
    end do
    call report_finalize()
end program test_zunbdb3
