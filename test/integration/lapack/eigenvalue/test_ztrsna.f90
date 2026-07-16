program test_ztrsna
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrsna
    use ref_quad_lapack, only: zgehrd, zhseqr, ztrevc, ztrsna
    implicit none
    integer, parameter :: ns(*) = [10, 20]
    integer :: i, n, info, lwork, m_real, m_r, m_g, j
    complex(ep), allocatable :: A(:,:), tau(:), work(:), T(:,:), Z(:,:), VL(:,:), VR(:,:), W(:)
    real(ep),    allocatable :: rwork(:), S_r(:), SEP_r(:), S_g(:), SEP_g(:)
    logical,     allocatable :: sel(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label
    call report_init('ztrsna', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 108201 + 47 * i)
        allocate(tau(n-1), T(n, n), Z(1, 1), W(n), sel(n))
        call zgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        T = A
        call zhseqr('S', 'N', n, 1, n, T, n, W, Z, 1, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhseqr('S', 'N', n, 1, n, T, n, W, Z, 1, work, lwork, info)
        deallocate(work)
        allocate(VL(n, n), VR(n, n), work(2*n), rwork(n))
        VL = (0.0_ep, 0.0_ep); do j = 1, n; VL(j, j) = (1.0_ep, 0.0_ep); end do
        VR = VL
        sel = .true.
        call ztrevc('B', 'A', sel, n, T, n, VL, n, VR, n, n, m_real, work, rwork, info)
        deallocate(work, rwork)
        allocate(S_r(n), SEP_r(n), S_g(n), SEP_g(n))
        allocate(work(n*(n+1)), rwork(n))
        call ztrsna('B', 'A', sel, n, T, n, VL, n, VR, n, S_r, SEP_r, n, m_r, &
                    work, n, rwork, info)
        deallocate(work, rwork)
        call target_ztrsna('B', 'A', sel, n, T, n, VL, n, VR, n, S_g, SEP_g, n, m_g, info)
        err = max(max_rel_err_vec(S_g(1:m_g), S_r(1:m_r)), max_rel_err_vec(SEP_g(1:m_g), SEP_r(1:m_r)))
        tol = 100.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(tau, T, Z, W, sel, VL, VR, S_r, SEP_r, S_g, SEP_g)
    end do
    call report_finalize()
end program test_ztrsna
