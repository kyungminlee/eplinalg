! dtrevc computes left and right eigenvectors of a quasi-triangular
! Schur matrix T from dgehrd + dhseqr. HOWMNY='A' (all), SIDE='R'.
program test_dtrevc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtrevc
    use ref_quad_lapack, only: dgehrd, dhseqr, dtrevc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_ref, m_got
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: T(:,:), WR(:), WI(:), Z(:,:)
    real(ep), allocatable :: VL(:,:), VR_ref(:,:), VR_got(:,:)
    logical,  allocatable :: sel_ref(:), sel_got(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtrevc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 350001 + 47 * i)
        allocate(tau(n-1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(T(n,n), WR(n), WI(n), Z(1,1))
        T = A
        call dhseqr('E', 'N', n, 1, n, T, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('E', 'N', n, 1, n, T, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        allocate(VL(1, n), VR_ref(n, n), VR_got(n, n))
        allocate(sel_ref(n), sel_got(n))
        sel_ref = .false.; sel_got = .false.
        allocate(work(3*n))
        call dtrevc('R', 'A', sel_ref, n, T, n, VL, 1, VR_ref, n, n, m_ref, work, info)
        call target_dtrevc('R', 'A', sel_got, n, T, n, VL, 1, VR_got, n, n, m_got, info)
        err = max_rel_err_mat(VR_got, VR_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, T, WR, WI, Z, VL, VR_ref, VR_got, sel_ref, sel_got, work)
    end do
    call report_finalize()
end program test_dtrevc
