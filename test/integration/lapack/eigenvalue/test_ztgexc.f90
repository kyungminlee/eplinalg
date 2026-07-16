! ztgexc: complex generalized Schur exchange. Smoke test using a
! proper generalized Schur form via QR(B) + zunmqr(A) + zgghrd + zhgeqz
! (mirrors test_dtgexc).
program test_ztgexc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgexc
    use ref_quad_lapack, only: zgeqrf, zunmqr, zgghrd, zhgeqz, ztgexc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, ifst, ilst, j
    complex(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    complex(ep), allocatable :: A_r(:,:), B_r(:,:), Q_r(:,:), Z_r(:,:)
    complex(ep), allocatable :: A_g(:,:), B_g(:,:), Q_g(:,:), Z_g(:,:)
    complex(ep), allocatable :: alpha(:), beta(:)
    real(ep), allocatable :: rwork(:), m_r(:), m_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgexc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 141001 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 141011 + 47 * i)
        allocate(tau(n))
        call zgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(n, n, B, n, tau, work, lwork, info)
        call zunmqr('L', 'C', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = (0.0_ep, 0.0_ep); end do
        allocate(A_r(n,n), B_r(n,n), Q_r(n,n), Z_r(n,n))
        allocate(A_g(n,n), B_g(n,n), Q_g(n,n), Z_g(n,n))
        allocate(alpha(n), beta(n), rwork(n), m_r(n), m_g(n))
        A_r = A; B_r = B
        Q_r = (0.0_ep, 0.0_ep); Z_r = (0.0_ep, 0.0_ep)
        do j = 1, n; Q_r(j,j) = (1.0_ep, 0.0_ep); Z_r(j,j) = (1.0_ep, 0.0_ep); end do
        call zgghrd('I', 'I', n, 1, n, A_r, n, B_r, n, Q_r, n, Z_r, n, info)
        call zhgeqz('S', 'V', 'V', n, 1, n, A_r, n, B_r, n, alpha, beta, &
                    Q_r, n, Z_r, n, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zhgeqz('S', 'V', 'V', n, 1, n, A_r, n, B_r, n, alpha, beta, &
                    Q_r, n, Z_r, n, work, lwork, rwork, info)
        deallocate(work)
        A_g = A_r; B_g = B_r; Q_g = Q_r; Z_g = Z_r
        ifst = n; ilst = 1
        call ztgexc(.true., .true., n, A_r, n, B_r, n, Q_r, n, Z_r, n, &
                    ifst, ilst, info)
        ifst = n; ilst = 1
        call target_ztgexc(.true., .true., n, A_g, n, B_g, n, Q_g, n, Z_g, n, &
                           ifst, ilst, info)
        do j = 1, n
            m_r(j) = merge(abs(A_r(j,j)) / abs(B_r(j,j)), huge(1.0_ep), &
                           abs(B_r(j,j)) > tiny(1.0_ep))
            m_g(j) = merge(abs(A_g(j,j)) / abs(B_g(j,j)), huge(1.0_ep), &
                           abs(B_g(j,j)) > tiny(1.0_ep))
        end do
        call sort_desc(m_r, n); call sort_desc(m_g, n)
        err = max_rel_err_vec(m_g, m_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, A_r, B_r, Q_r, Z_r, A_g, B_g, Q_g, Z_g, &
                   alpha, beta, rwork, m_r, m_g)
    end do
    call report_finalize()
contains
    subroutine sort_desc(x, m)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: m
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, m - 1
            do jj = ii + 1, m
                if (x(ii) < x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_ztgexc
