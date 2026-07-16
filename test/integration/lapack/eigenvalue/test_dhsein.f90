! dhsein: inverse iteration eigenvectors of an upper Hessenberg matrix
! given approximate eigenvalues. SIDE='R', EIGSRC='N', INITV='N'.
! Inverse-iteration eigenvectors are not unique (sign and, at lower
! precision, the iteration may converge to a different vector in the
! same invariant subspace), so we check the eigenvector residual
! ||H*v - lambda*v|| / (||H||*||v||) instead of comparing entry-wise.
program test_dhsein
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dhsein
    use ref_quad_lapack, only: dgehrd, dhseqr
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, m_g, j, ii
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: H(:,:), WR(:), WI(:), Z(:,:)
    real(ep), allocatable :: VL(:,:), VR(:,:), Hfull(:,:)
    real(ep), allocatable :: u(:), wv(:), Hu(:), Hw(:)
    integer,  allocatable :: ifaill(:), ifailr(:)
    logical,  allocatable :: select(:)
    real(ep) :: wopt(1), err, tol, hnorm, vnorm, rnorm
    character(len=48) :: label

    call report_init('dhsein', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 112001 + 47 * i)
        allocate(tau(n - 1))
        call dgehrd(n, 1, n, A, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dgehrd(n, 1, n, A, n, tau, work, lwork, info)
        deallocate(work)
        allocate(H(n,n), WR(n), WI(n), Z(1,1))
        H = A
        call dhseqr('E', 'N', n, 1, n, H, n, WR, WI, Z, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dhseqr('E', 'N', n, 1, n, H, n, WR, WI, Z, 1, work, lwork, info)
        deallocate(work)
        ! Restore H to its Hessenberg form (overwritten by dhseqr).
        H = A
        allocate(VL(1, n), VR(n, n), ifaill(n), ifailr(n), select(n))
        select = .true.
        VR = 0.0_ep
        call target_dhsein('R', 'N', 'N', select, n, H, n, WR, WI, &
                           VL, 1, VR, n, n, m_g, ifaill, ifailr, info)
        ! Reconstruct the full Hessenberg matrix (only upper triangle +
        ! first subdiagonal are used; lower trapezoid in H holds reflectors).
        allocate(Hfull(n, n), u(n), wv(n), Hu(n), Hw(n))
        Hfull = 0.0_ep
        do j = 1, n
            do ii = 1, min(j + 1, n)
                Hfull(ii, j) = H(ii, j)
            end do
        end do
        hnorm = maxval(sum(abs(Hfull), dim = 1))
        err = 0.0_ep
        j = 1
        do while (j <= n)
            if (j < n .and. WI(j) /= 0.0_ep) then
                u = VR(:, j); wv = VR(:, j + 1)
                Hu = matmul(Hfull, u); Hw = matmul(Hfull, wv)
                vnorm = sqrt(sum(u * u) + sum(wv * wv))
                rnorm = sqrt(sum((Hu - WR(j) * u + WI(j) * wv)**2) &
                          +  sum((Hw - WI(j) * u - WR(j) * wv)**2))
                err = max(err, rnorm / max(hnorm * vnorm, tiny(1.0_ep)))
                j = j + 2
            else
                u = VR(:, j)
                Hu = matmul(Hfull, u)
                vnorm = sqrt(sum(u * u))
                rnorm = sqrt(sum((Hu - WR(j) * u)**2))
                err = max(err, rnorm / max(hnorm * vnorm, tiny(1.0_ep)))
                j = j + 1
            end if
        end do
        deallocate(Hfull, u, wv, Hu, Hw)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, tau, H, WR, WI, Z, VL, VR, ifaill, ifailr, select)
    end do
    call report_finalize()
end program test_dhsein
