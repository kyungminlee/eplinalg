! dstemr: MRRR tridiagonal eigensolver with try-rac.
program test_dstemr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use target_lapack,   only: target_name, target_eps, target_dstemr
    use ref_quad_lapack, only: dstemr
    implicit none

    integer, parameter :: ns(*) = [10, 32, 64]
    integer :: i, n, m_r, m_g, info, lwork, liwork, j, nzc
    real(ep), allocatable :: D0(:), E0(:), D_r(:), D_g(:), E_r(:), E_g(:)
    real(ep), allocatable :: W_r(:), W_g(:), Z(:,:), work(:)
    integer, allocatable :: isuppz(:), iwork(:)
    real(ep) :: wopt(1), err, tol
    logical :: tryrac_r, tryrac_g
    integer :: iopt(1)
    character(len=48) :: label

    call report_init('dstemr', target_name)
    do i = 1, size(ns)
        n = ns(i); nzc = 0
        allocate(D0(n), E0(max(1, n)), D_r(n), D_g(n), E_r(max(1, n)), E_g(max(1, n)))
        allocate(W_r(n), W_g(n), Z(1, 1), isuppz(2*n))
        do j = 1, n; D0(j) = real(j, ep); end do
        do j = 1, n-1; E0(j) = 0.5_ep; end do
        E0(max(1, n)) = 0.0_ep
        D_r = D0; D_g = D0; E_r = E0; E_g = E0
        tryrac_r = .true.; tryrac_g = .true.
        call dstemr('N', 'A', n, D_r, E_r, 0.0_ep, 0.0_ep, 0, 0, &
                    m_r, W_r, Z, 1, nzc, isuppz, tryrac_r, wopt, -1, iopt, -1, info)
        lwork = max(1, int(wopt(1))); liwork = max(1, iopt(1))
        allocate(work(lwork), iwork(liwork))
        call dstemr('N', 'A', n, D_r, E_r, 0.0_ep, 0.0_ep, 0, 0, &
                    m_r, W_r, Z, 1, nzc, isuppz, tryrac_r, work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        call target_dstemr('N', 'A', n, D_g, E_g, 0.0_ep, 0.0_ep, 0, 0, &
                           m_g, W_g, Z, 1, nzc, isuppz, tryrac_g, info)
        err = max_rel_err_vec(W_g(1:m_g), W_r(1:m_r))
        if (m_r /= m_g) err = max(err, 1.0_ep)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(D0, E0, D_r, D_g, E_r, E_g, W_r, W_g, Z, isuppz)
    end do
    call report_finalize()
end program test_dstemr
