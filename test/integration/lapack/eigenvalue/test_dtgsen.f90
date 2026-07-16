! dtgsen: reorder generalized Schur (A,B). Smoke test using upper-tri
! pair (no swap candidates beyond adjacent rows). Compare sorted moduli
! after selecting the first half of eigenvalues.
program test_dtgsen
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtgsen
    use ref_quad_lapack, only: dgeqrf, dormqr, dgghrd, dhgeqz, dtgsen
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j, mout
    real(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    real(ep), allocatable :: A_r(:,:), B_r(:,:), Q_r(:,:), Z_r(:,:)
    real(ep), allocatable :: A_g(:,:), B_g(:,:), Q_g(:,:), Z_g(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:)
    real(ep), allocatable :: ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: ar(:), ai(:), be(:), m_r(:), m_g(:)
    logical, allocatable :: select(:)
    integer, allocatable :: iwork(:)
    real(ep) :: wopt(1), pl_r, pr_r, dif_r(2), pl_g, pr_g, dif_g(2)
    integer :: iwopt(1), liwork
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dtgsen', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 144001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 144011 + 47 * i)
        allocate(tau(n))
        call dgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(n, n, B, n, tau, work, lwork, info)
        call dormqr('L', 'T', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = 0.0_ep; end do
        allocate(A_r(n,n), B_r(n,n), Q_r(n,n), Z_r(n,n))
        allocate(A_g(n,n), B_g(n,n), Q_g(n,n), Z_g(n,n))
        allocate(ar(n), ai(n), be(n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(m_r(n), m_g(n), select(n))
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
        select = .false.
        do j = 1, n / 2; select(j) = .true.; end do
        call dtgsen(0, .true., .true., select, n, A_r, n, B_r, n, &
                    ar_r, ai_r, be_r, Q_r, n, Z_r, n, mout, &
                    pl_r, pr_r, dif_r, wopt, -1, iwopt, -1, info)
        lwork = max(1, int(wopt(1))); liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dtgsen(0, .true., .true., select, n, A_r, n, B_r, n, &
                    ar_r, ai_r, be_r, Q_r, n, Z_r, n, mout, &
                    pl_r, pr_r, dif_r, work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        call target_dtgsen(0, .true., .true., select, n, A_g, n, B_g, n, &
                           ar_g, ai_g, be_g, Q_g, n, Z_g, n, mout, &
                           pl_g, pr_g, dif_g, info)
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
                   ar, ai, be, ar_r, ai_r, be_r, ar_g, ai_g, be_g, &
                   m_r, m_g, select)
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
end program test_dtgsen
