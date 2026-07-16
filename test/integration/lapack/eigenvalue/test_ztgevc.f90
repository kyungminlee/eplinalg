! ztgevc: complex generalized Schur eigenvectors. SIDE='R', HOWMNY='A'.
program test_ztgevc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgevc
    use ref_quad_lapack, only: zgeqrf, zunmqr, zgghrd, zhgeqz, ztgevc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_r, m_g, j
    complex(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    complex(ep), allocatable :: S(:,:), P(:,:), Q(:,:), Z(:,:), a_(:), be(:)
    complex(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    real(ep),    allocatable :: rwork(:)
    logical,     allocatable :: select(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgevc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 133001 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 133011 + 47 * i)
        allocate(tau(n))
        call zgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(n, n, B, n, tau, work, lwork, info)
        call zunmqr('L', 'C', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = (0.0_ep, 0.0_ep); end do
        allocate(S(n,n), P(n,n), Q(n,n), Z(n,n), a_(n), be(n), rwork(n))
        S = A; P = B; Q = (0.0_ep, 0.0_ep); Z = (0.0_ep, 0.0_ep)
        do j = 1, n
            Q(j,j) = (1.0_ep, 0.0_ep); Z(j,j) = (1.0_ep, 0.0_ep)
        end do
        call zgghrd('N', 'N', n, 1, n, S, n, P, n, Q, n, Z, n, info)
        call zhgeqz('S', 'N', 'N', n, 1, n, S, n, P, n, a_, be, &
                    Q, n, Z, n, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zhgeqz('S', 'N', 'N', n, 1, n, S, n, P, n, a_, be, &
                    Q, n, Z, n, work, lwork, rwork, info)
        deallocate(work)
        allocate(VL(1,n), VR_r(n,n), VR_g(n,n), select(n))
        allocate(work(2*n))
        if (size(rwork) < 2*n) then
            deallocate(rwork); allocate(rwork(2*n))
        end if
        select = .false.
        call ztgevc('R', 'A', select, n, S, n, P, n, VL, 1, VR_r, n, n, &
                    m_r, work, rwork, info)
        deallocate(work)
        call target_ztgevc('R', 'A', select, n, S, n, P, n, VL, 1, VR_g, n, n, &
                           m_g, info)
        err = max_rel_err_mat_z(VR_g, VR_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, S, P, Q, Z, a_, be, rwork, VL, VR_r, VR_g, select)
    end do
    call report_finalize()
end program test_ztgevc
