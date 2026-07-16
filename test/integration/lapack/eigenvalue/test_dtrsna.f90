program test_dtrsna
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrsna
    use ref_quad_lapack, only: dgehrd, dhseqr, dtrevc, dtrsna
    implicit none
    integer, parameter :: ns(*) = [10, 20]
    integer :: i, n, info, lwork, m_real, m_r, m_g
    real(ep), allocatable :: A(:,:), tau(:), work(:), T(:,:), Z(:,:)
    real(ep), allocatable :: WR(:), WI(:), VL(:,:), VR(:,:)
    real(ep), allocatable :: S_r(:), SEP_r(:), S_g(:), SEP_g(:)
    integer,  allocatable :: iwork(:)
    logical,  allocatable :: sel(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label
    call report_init('dtrsna', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 108101 + 47 * i)
        allocate(tau(n-1), T(n, n), Z(1, 1), WR(n), WI(n), sel(n))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        T = A
        call dhseqr('S', 'N', n, 1, n, T, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('S', 'N', n, 1, n, T, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        allocate(VL(n, n), VR(n, n), work(3*n))
        VL = 0.0_ep; do m_r = 1, n; VL(m_r, m_r) = 1.0_ep; end do
        VR = VL
        sel = .true.
        call dtrevc('B', 'A', sel, n, T, n, VL, n, VR, n, n, m_real, work, info)
        deallocate(work)
        allocate(S_r(n), SEP_r(n), S_g(n), SEP_g(n))
        allocate(work(n*(n+6)), iwork(2*(n-1)))
        call dtrsna('B', 'A', sel, n, T, n, VL, n, VR, n, S_r, SEP_r, n, m_r, &
                    work, n, iwork, info)
        deallocate(work, iwork)
        call target_dtrsna('B', 'A', sel, n, T, n, VL, n, VR, n, S_g, SEP_g, n, m_g, info)
        err = max(max_rel_err_vec(S_g(1:m_g), S_r(1:m_r)), max_rel_err_vec(SEP_g(1:m_g), SEP_r(1:m_r)))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(tau, T, Z, WR, WI, sel, VL, VR, S_r, SEP_r, S_g, SEP_g)
    end do
    call report_finalize()
end program test_dtrsna
