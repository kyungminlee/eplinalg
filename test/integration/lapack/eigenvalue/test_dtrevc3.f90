! dtrevc3: blocked TREVC for real Schur eigenvectors. SIDE='R',
! HOWMNY='A'. Right eigenvectors of T computed.
program test_dtrevc3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrevc3
    use ref_quad_lapack, only: dgehrd, dhseqr, dtrevc3
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_r, m_g
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: T(:,:), WR(:), WI(:), Z(:,:)
    real(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    logical,  allocatable :: select(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtrevc3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 110001 + 47 * i)
        allocate(tau(n - 1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(T(n,n), WR(n), WI(n), Z(1,1))
        T = A
        call dhseqr('S', 'N', n, 1, n, T, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('S', 'N', n, 1, n, T, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        allocate(VL(1, n), VR_r(n, n), VR_g(n, n), select(n))
        select = .false.
        call dtrevc3('R', 'A', select, n, T, n, VL, 1, VR_r, n, n, m_r, &
                     wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dtrevc3('R', 'A', select, n, T, n, VL, 1, VR_r, n, n, m_r, &
                     work, lwork, info)
        deallocate(work)
        call target_dtrevc3('R', 'A', select, n, T, n, VL, 1, VR_g, n, n, &
                            m_g, info)
        err = max_rel_err_mat(VR_g, VR_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, T, WR, WI, Z, VL, VR_r, VR_g, select)
    end do
    call report_finalize()
end program test_dtrevc3
