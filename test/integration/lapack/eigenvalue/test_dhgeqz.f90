! dhgeqz: QZ algorithm for generalized eig of (H,T) Hess-tri pair.
! Setup: dgghrd reduces (A,B) to (H,T); dhgeqz computes eigvals.
program test_dhgeqz
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dhgeqz
    use ref_quad_lapack, only: dgeqrf, dormqr, dgghrd, dhgeqz
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j
    real(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    real(ep), allocatable :: H_r(:,:), T_r(:,:), H_g(:,:), T_g(:,:), Q(:,:), Z(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:), ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: m_r(:), m_g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dhgeqz', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 130001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 130011 + 47 * i)
        ! QR factorize B then apply Q^T to A so B becomes upper triangular.
        allocate(tau(n))
        call dgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(n, n, B, n, tau, work, lwork, info)
        call dormqr('L', 'T', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        ! Zero strict-lower triangle of B.
        do j = 1, n
            B(j+1:n, j) = 0.0_ep
        end do
        ! Reduce A to Hessenberg, B stays triangular.
        allocate(H_r(n,n), T_r(n,n), H_g(n,n), T_g(n,n), Q(n,n), Z(n,n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(m_r(n), m_g(n))
        H_r = A; T_r = B; H_g = A; T_g = B
        Q = 0.0_ep; Z = 0.0_ep
        do j = 1, n; Q(j,j) = 1.0_ep; Z(j,j) = 1.0_ep; end do
        call dgghrd('N', 'N', n, 1, n, H_r, n, T_r, n, Q, n, Z, n, info)
        H_g = H_r; T_g = T_r
        ! Compute eigvals via QZ.
        call dhgeqz('E', 'N', 'N', n, 1, n, H_r, n, T_r, n, ar_r, ai_r, be_r, &
                    Q, n, Z, n, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dhgeqz('E', 'N', 'N', n, 1, n, H_r, n, T_r, n, ar_r, ai_r, be_r, &
                    Q, n, Z, n, work, lwork, info)
        deallocate(work)
        call target_dhgeqz('E', 'N', 'N', n, 1, n, H_g, n, T_g, n, ar_g, ai_g, be_g, &
                           Q, n, Z, n, info)
        do j = 1, n
            m_r(j) = merge(sqrt(ar_r(j)**2 + ai_r(j)**2) / abs(be_r(j)), &
                           huge(1.0_ep), abs(be_r(j)) > tiny(1.0_ep))
            m_g(j) = merge(sqrt(ar_g(j)**2 + ai_g(j)**2) / abs(be_g(j)), &
                           huge(1.0_ep), abs(be_g(j)) > tiny(1.0_ep))
        end do
        call sort_desc(m_r, n); call sort_desc(m_g, n)
        err = max_rel_err_vec(m_g, m_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, H_r, T_r, H_g, T_g, Q, Z, &
                   ar_r, ai_r, be_r, ar_g, ai_g, be_g, m_r, m_g)
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
end program test_dhgeqz
