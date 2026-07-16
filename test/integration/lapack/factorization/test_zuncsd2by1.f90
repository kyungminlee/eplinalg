! zuncsd2by1: 2-by-1 CS decomposition of [X11; X21] (complex unitary block).
program test_zuncsd2by1
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zuncsd2by1
    use ref_quad_lapack, only: zgeqrf, zungqr, zuncsd2by1
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, p, q, info, lwork, lrwork
    complex(ep), allocatable :: A(:,:), tau(:), work(:)
    complex(ep), allocatable :: X11r(:,:), X21r(:,:), X11g(:,:), X21g(:,:)
    complex(ep), allocatable :: U1(:,:), U2(:,:), V1t(:,:)
    real(ep), allocatable :: theta_r(:), theta_g(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: ropt(1), err, tol
    integer, allocatable :: iwork(:)
    character(len=48) :: label

    call report_init('zuncsd2by1', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = m/2
        call gen_matrix_complex(m, m, A, seed = 24251 + 79 * i)
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
        allocate(U1(p, p), U2(m-p, m-p), V1t(q, q), &
                 theta_r(q), theta_g(q), iwork(m-q))
        call zuncsd2by1('Y', 'Y', 'Y', m, p, q, X11r, p, X21r, m-p, &
                        theta_r, U1, p, U2, m-p, V1t, q, &
                        wopt, -1, ropt, -1, iwork, info)
        lwork  = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        lrwork = max(1, int(ropt(1)));            allocate(rwork(lrwork))
        call zuncsd2by1('Y', 'Y', 'Y', m, p, q, X11r, p, X21r, m-p, &
                        theta_r, U1, p, U2, m-p, V1t, q, &
                        work, lwork, rwork, lrwork, iwork, info)
        deallocate(work, rwork)
        call target_zuncsd2by1('Y', 'Y', 'Y', m, p, q, X11g, p, X21g, m-p, &
                               theta_g, U1, p, U2, m-p, V1t, q, info)
        err = max_rel_err_vec(theta_g, theta_r)
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X21r, X11g, X21g, U1, U2, V1t, theta_r, theta_g, iwork)
    end do
    call report_finalize()
end program test_zuncsd2by1
