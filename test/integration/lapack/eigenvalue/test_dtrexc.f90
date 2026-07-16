! dtrexc: exchange diagonal blocks of a Schur form. Reduce A to Schur,
! swap blocks, compare; the eigenvalue spectrum must be invariant.
program test_dtrexc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrexc
    use ref_quad_lapack, only: dgehrd, dhseqr, dtrexc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, ifst, ilst, j
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: T_r(:,:), T_g(:,:), Q_r(:,:), Q_g(:,:)
    real(ep), allocatable :: WR(:), WI(:), Z(:,:)
    real(ep), allocatable :: diag_r(:), diag_g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtrexc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 106001 + 47 * i)
        allocate(tau(n - 1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(T_r(n, n), T_g(n, n), Q_r(n, n), Q_g(n, n), WR(n), WI(n), Z(1,1))
        T_r = A; T_g = A
        call dhseqr('S', 'N', n, 1, n, T_r, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('S', 'N', n, 1, n, T_r, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        T_g = T_r
        Q_r = 0.0_ep; Q_g = 0.0_ep
        do j = 1, n; Q_r(j, j) = 1.0_ep; Q_g(j, j) = 1.0_ep; end do
        ifst = n; ilst = 1
        allocate(work(n))
        call dtrexc('V', n, T_r, n, Q_r, n, ifst, ilst, work, info)
        deallocate(work)
        ifst = n; ilst = 1
        call target_dtrexc('V', n, T_g, n, Q_g, n, ifst, ilst, info)
        ! The diagonal of T should hold the same eigenvalues (modulo
        ! 2x2 blocks for complex pairs); compare sorted diagonals.
        allocate(diag_r(n), diag_g(n))
        do j = 1, n
            diag_r(j) = T_r(j, j)
            diag_g(j) = T_g(j, j)
        end do
        call sort_asc(diag_r, n); call sort_asc(diag_g, n)
        err = max_rel_err_vec(diag_g, diag_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, T_r, T_g, Q_r, Q_g, WR, WI, Z, diag_r, diag_g)
    end do
    call report_finalize()
contains
    subroutine sort_asc(x, m)
        real(ep), intent(inout) :: x(:)
        integer,  intent(in)    :: m
        integer :: ii, jj
        real(ep) :: tt
        do ii = 1, m - 1
            do jj = ii + 1, m
                if (x(ii) > x(jj)) then
                    tt = x(ii); x(ii) = x(jj); x(jj) = tt
                end if
            end do
        end do
    end subroutine
end program test_dtrexc
