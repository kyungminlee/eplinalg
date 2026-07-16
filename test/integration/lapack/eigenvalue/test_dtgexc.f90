! dtgexc: exchange diagonal blocks of generalized Schur (A,B). Smoke
! test — set up an upper-triangular pair (no swap candidates beyond
! adjacent rows), call routine, check INFO and a sorted invariant.
program test_dtgexc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtgexc
    use ref_quad_lapack, only: dgeqrf, dormqr, dgghrd, dhgeqz, dtgexc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, ifst, ilst, j
    real(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    real(ep), allocatable :: A_r(:,:), B_r(:,:), Q_r(:,:), Z_r(:,:)
    real(ep), allocatable :: A_g(:,:), B_g(:,:), Q_g(:,:), Z_g(:,:)
    real(ep), allocatable :: ar(:), ai(:), be(:), m_r(:), m_g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtgexc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 140001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 140011 + 47 * i)
        allocate(tau(n))
        call dgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(n, n, B, n, tau, work, lwork, info)
        call dormqr('L', 'T', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = 0.0_ep; end do
        allocate(A_r(n,n), B_r(n,n), Q_r(n,n), Z_r(n,n))
        allocate(A_g(n,n), B_g(n,n), Q_g(n,n), Z_g(n,n))
        allocate(ar(n), ai(n), be(n), m_r(n), m_g(n))
        A_r = A; B_r = B; Q_r = 0.0_ep; Z_r = 0.0_ep
        do j = 1, n; Q_r(j,j) = 1.0_ep; Z_r(j,j) = 1.0_ep; end do
        call dgghrd('I', 'I', n, 1, n, A_r, n, B_r, n, Q_r, n, Z_r, n, info)
        call dhgeqz('S', 'V', 'V', n, 1, n, A_r, n, B_r, n, ar, ai, be, &
                    Q_r, n, Z_r, n, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dhgeqz('S', 'V', 'V', n, 1, n, A_r, n, B_r, n, ar, ai, be, &
                    Q_r, n, Z_r, n, work, lwork, info)
        deallocate(work)
        A_g = A_r; B_g = B_r; Q_g = Q_r; Z_g = Z_r
        ifst = n; ilst = 1
        allocate(work(4*n+16))
        call dtgexc(.true., .true., n, A_r, n, B_r, n, Q_r, n, Z_r, n, &
                    ifst, ilst, work, 4*n+16, info)
        deallocate(work)
        ifst = n; ilst = 1
        call target_dtgexc(.true., .true., n, A_g, n, B_g, n, Q_g, n, Z_g, n, &
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
                   ar, ai, be, m_r, m_g)
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
end program test_dtgexc
