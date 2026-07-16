! zunbdb: complex variant of dorbdb.
program test_zunbdb
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zunbdb
    use ref_quad_lapack, only: zunbdb
    implicit none

    integer, parameter :: ms(*) = [12, 24]
    integer :: i, m, p, q, info, lwork
    complex(ep), allocatable :: X11r(:,:), X12r(:,:), X21r(:,:), X22r(:,:)
    complex(ep), allocatable :: X11g(:,:), X12g(:,:), X21g(:,:), X22g(:,:)
    real(ep),    allocatable :: theta_r(:), theta_g(:), phi_r(:), phi_g(:)
    complex(ep), allocatable :: tp1r(:), tp2r(:), tq1r(:), tq2r(:), work(:)
    complex(ep), allocatable :: tp1g(:), tp2g(:), tq1g(:), tq2g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zunbdb', target_name)
    do i = 1, size(ms)
        m = ms(i); p = m/2; q = m/2
        allocate(X11r(p, q), X12r(p, m-q), X21r(m-p, q), X22r(m-p, m-q))
        allocate(X11g(p, q), X12g(p, m-q), X21g(m-p, q), X22g(m-p, m-q))
        allocate(theta_r(q), theta_g(q), phi_r(max(1, q-1)), phi_g(max(1, q-1)))
        allocate(tp1r(p), tp2r(m-p), tq1r(q), tq2r(m-q))
        allocate(tp1g(p), tp2g(m-p), tq1g(q), tq2g(m-q))
        block
            complex(ep), allocatable :: A(:,:)
            integer :: ii, jj
            call gen_matrix_complex(m, m, A, seed = 24051 + 71 * i)
            do ii = 1, p
                do jj = 1, q;   X11r(ii, jj) = 0.1_ep * A(ii, jj); end do
                do jj = 1, m-q; X12r(ii, jj) = 0.1_ep * A(ii, q+jj); end do
            end do
            do ii = 1, m-p
                do jj = 1, q;   X21r(ii, jj) = 0.1_ep * A(p+ii, jj); end do
                do jj = 1, m-q; X22r(ii, jj) = 0.1_ep * A(p+ii, q+jj); end do
            end do
            do ii = 1, p; X11r(ii, ii) = X11r(ii, ii) + (1.0_ep, 0.0_ep); end do
            do ii = 1, m-p; X22r(ii, ii) = X22r(ii, ii) + (1.0_ep, 0.0_ep); end do
            deallocate(A)
        end block
        X11g = X11r; X12g = X12r; X21g = X21r; X22g = X22r
        call zunbdb('N', 'O', m, p, q, X11r, p, X12r, p, X21r, m-p, X22r, m-p, &
                    theta_r, phi_r, tp1r, tp2r, tq1r, tq2r, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zunbdb('N', 'O', m, p, q, X11r, p, X12r, p, X21r, m-p, X22r, m-p, &
                    theta_r, phi_r, tp1r, tp2r, tq1r, tq2r, work, lwork, info)
        deallocate(work)
        call target_zunbdb('N', 'O', m, p, q, X11g, p, X12g, p, X21g, m-p, X22g, m-p, &
                           theta_g, phi_g, tp1g, tp2g, tq1g, tq2g, info)
        err = max_rel_err_vec(theta_g, theta_r)
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0)') 'm=', m
        call report_case(trim(label), err, tol)
        deallocate(X11r, X12r, X21r, X22r, X11g, X12g, X21g, X22g, &
                   theta_r, theta_g, phi_r, phi_g, &
                   tp1r, tp2r, tq1r, tq2r, tp1g, tp2g, tq1g, tq2g)
    end do
    call report_finalize()
end program test_zunbdb
