! dtrsen: condition number / span for selected eigenvalues. JOB='N',
! COMPQ='N' — no eigenvector accumulation; just reorder eigenvalues
! marked by SELECT to the top of T.
program test_dtrsen
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrsen
    use ref_quad_lapack, only: dgehrd, dhseqr, dtrsen
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, liwork, m_r, m_g, j, iwopt(1)
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: T_r(:,:), T_g(:,:), Q(:,:)
    real(ep), allocatable :: WR_r(:), WI_r(:), WR_g(:), WI_g(:), Z(:,:)
    integer,  allocatable :: iwork(:)
    logical,  allocatable :: select(:)
    real(ep), allocatable :: mods_r(:), mods_g(:)
    real(ep) :: wopt(1), s_r, sep_r, s_g, sep_g, err, tol
    character(len=48) :: label

    call report_init('dtrsen', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 108001 + 47 * i)
        allocate(tau(n - 1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(T_r(n, n), T_g(n, n), Q(1,1), WR_r(n), WI_r(n))
        allocate(WR_g(n), WI_g(n), Z(1,1), select(n))
        T_r = A
        call dhseqr('S', 'N', n, 1, n, T_r, n, WR_r, WI_r, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('S', 'N', n, 1, n, T_r, n, WR_r, WI_r, Z, 1, work, lwork, info)
        deallocate(work)
        T_g = T_r
        select = .false.
        do j = 1, n / 2
            select(j) = .true.
        end do
        call dtrsen('N', 'N', select, n, T_r, n, Q, 1, WR_r, WI_r, m_r, &
                    s_r, sep_r, wopt, -1, iwopt, -1, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dtrsen('N', 'N', select, n, T_r, n, Q, 1, WR_r, WI_r, m_r, &
                    s_r, sep_r, work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        call target_dtrsen('N', 'N', select, n, T_g, n, Q, 1, WR_g, WI_g, &
                           m_g, s_g, sep_g, info)
        allocate(mods_r(n), mods_g(n))
        do j = 1, n
            mods_r(j) = sqrt(WR_r(j)**2 + WI_r(j)**2)
            mods_g(j) = sqrt(WR_g(j)**2 + WI_g(j)**2)
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, T_r, T_g, Q, WR_r, WI_r, WR_g, WI_g, Z, &
                   select, mods_r, mods_g)
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
end program test_dtrsen
