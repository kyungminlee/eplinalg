! dtgevc: eigenvectors of generalized Schur (S,P). SIDE='R', HOWMNY='A'.
program test_dtgevc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtgevc
    use ref_quad_lapack, only: dgeqrf, dormqr, dgghrd, dhgeqz, dtgevc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_r, m_g, j
    real(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    real(ep), allocatable :: S(:,:), P(:,:), Q(:,:), Z(:,:)
    real(ep), allocatable :: ar(:), ai(:), be(:)
    real(ep), allocatable :: VL(:,:), VR_r(:,:), VR_g(:,:)
    logical,  allocatable :: select(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtgevc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 132001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 132011 + 47 * i)
        allocate(tau(n))
        call dgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(n, n, B, n, tau, work, lwork, info)
        call dormqr('L', 'T', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = 0.0_ep; end do
        allocate(S(n,n), P(n,n), Q(n,n), Z(n,n), ar(n), ai(n), be(n))
        S = A; P = B; Q = 0.0_ep; Z = 0.0_ep
        do j = 1, n; Q(j,j) = 1.0_ep; Z(j,j) = 1.0_ep; end do
        call dgghrd('N', 'N', n, 1, n, S, n, P, n, Q, n, Z, n, info)
        call dhgeqz('S', 'N', 'N', n, 1, n, S, n, P, n, ar, ai, be, &
                    Q, n, Z, n, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dhgeqz('S', 'N', 'N', n, 1, n, S, n, P, n, ar, ai, be, &
                    Q, n, Z, n, work, lwork, info)
        deallocate(work)
        allocate(VL(1,n), VR_r(n,n), VR_g(n,n), select(n), work(6*n))
        select = .false.
        call dtgevc('R', 'A', select, n, S, n, P, n, VL, 1, VR_r, n, n, &
                    m_r, work, info)
        deallocate(work)
        call target_dtgevc('R', 'A', select, n, S, n, P, n, VL, 1, VR_g, n, n, &
                           m_g, info)
        err = max_rel_err_mat(VR_g, VR_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, S, P, Q, Z, ar, ai, be, VL, VR_r, VR_g, select)
    end do
    call report_finalize()
end program test_dtgevc
