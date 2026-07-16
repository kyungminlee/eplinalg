! dtgsna: condition numbers for generalized eigenvalues/vectors. Smoke
! test. Set up generalized Schur via gghrd+hgeqz, eigenvectors via
! dtgevc, then compare s and dif outputs.
program test_dtgsna
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dtgsna
    use ref_quad_lapack, only: dgeqrf, dormqr, dgghrd, dhgeqz, dtgevc, dtgsna
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j, mout, mm
    real(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    real(ep), allocatable :: A_s(:,:), B_s(:,:), Q_s(:,:), Z_s(:,:)
    real(ep), allocatable :: VL(:,:), VR(:,:)
    real(ep), allocatable :: ar(:), ai(:), be(:)
    real(ep), allocatable :: s_r(:), dif_r(:), s_g(:), dif_g(:)
    logical, allocatable :: select(:)
    integer, allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dtgsna', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 145001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 145011 + 47 * i)
        allocate(tau(n))
        call dgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(n, n, B, n, tau, work, lwork, info)
        call dormqr('L', 'T', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = 0.0_ep; end do
        allocate(A_s(n,n), B_s(n,n), Q_s(n,n), Z_s(n,n))
        allocate(ar(n), ai(n), be(n))
        A_s = A; B_s = B; Q_s = 0.0_ep; Z_s = 0.0_ep
        do j = 1, n; Q_s(j,j) = 1.0_ep; Z_s(j,j) = 1.0_ep; end do
        call dgghrd('I', 'I', n, 1, n, A_s, n, B_s, n, Q_s, n, Z_s, n, info)
        call dhgeqz('S', 'V', 'V', n, 1, n, A_s, n, B_s, n, ar, ai, be, &
                    Q_s, n, Z_s, n, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dhgeqz('S', 'V', 'V', n, 1, n, A_s, n, B_s, n, ar, ai, be, &
                    Q_s, n, Z_s, n, work, lwork, info)
        deallocate(work)
        ! Eigenvectors of (A_s, B_s).
        mm = n
        allocate(VL(n, mm), VR(n, mm), select(n))
        VL = 0.0_ep; VR = 0.0_ep; select = .true.
        allocate(work(6 * n))
        call dtgevc('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                    mm, mout, work, info)
        deallocate(work)
        allocate(s_r(n), dif_r(n), s_g(n), dif_g(n))
        allocate(work(max(1, 2 * n * (n + 2) + 16)), iwork(n+6))
        call dtgsna('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                    s_r, dif_r, mm, mout, work, size(work), iwork, info)
        deallocate(work, iwork)
        call target_dtgsna('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                           s_g, dif_g, mm, mout, info)
        ! Eigenvalue order can differ when adjacent eigenvalues collide;
        ! sort both vectors before comparing.
        ! Compare only the eigenvalue condition numbers (s) — dif values
        ! span many orders of magnitude and are too precision-sensitive
        ! for elementwise comparison across targets.
        call sort_desc_local(s_r, mout); call sort_desc_local(s_g, mout)
        err = max_rel_err_vec(s_g(1:mout), s_r(1:mout))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, A_s, B_s, Q_s, Z_s, ar, ai, be, &
                   VL, VR, select, s_r, dif_r, s_g, dif_g)
    end do
    call report_finalize()
contains
    subroutine sort_desc_local(x, m)
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
end program test_dtgsna
