! ztgsna: complex generalized eigenvalue/vector condition numbers.
! Smoke test using a proper generalized Schur form via QR(B) +
! zunmqr(A) + zgghrd + zhgeqz, eigenvectors via ztgevc.
program test_ztgsna
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgsna
    use ref_quad_lapack, only: zgeqrf, zunmqr, zgghrd, zhgeqz, ztgevc, ztgsna
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j, mout, mm
    complex(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    complex(ep), allocatable :: A_s(:,:), B_s(:,:), Q_s(:,:), Z_s(:,:)
    complex(ep), allocatable :: VL(:,:), VR(:,:), alpha(:), beta(:)
    real(ep), allocatable :: s_r(:), dif_r(:), s_g(:), dif_g(:), rwork(:)
    logical, allocatable :: select(:)
    integer, allocatable :: iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgsna', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 145101 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 145111 + 47 * i)
        allocate(tau(n))
        call zgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(n, n, B, n, tau, work, lwork, info)
        call zunmqr('L', 'C', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = (0.0_ep, 0.0_ep); end do
        allocate(A_s(n,n), B_s(n,n), Q_s(n,n), Z_s(n,n))
        allocate(alpha(n), beta(n), rwork(2*n))
        A_s = A; B_s = B
        Q_s = (0.0_ep, 0.0_ep); Z_s = (0.0_ep, 0.0_ep)
        do j = 1, n; Q_s(j,j) = (1.0_ep, 0.0_ep); Z_s(j,j) = (1.0_ep, 0.0_ep); end do
        call zgghrd('I', 'I', n, 1, n, A_s, n, B_s, n, Q_s, n, Z_s, n, info)
        call zhgeqz('S', 'V', 'V', n, 1, n, A_s, n, B_s, n, alpha, beta, &
                    Q_s, n, Z_s, n, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zhgeqz('S', 'V', 'V', n, 1, n, A_s, n, B_s, n, alpha, beta, &
                    Q_s, n, Z_s, n, work, lwork, rwork, info)
        deallocate(work)
        mm = n
        allocate(VL(n, mm), VR(n, mm), select(n))
        VL = (0.0_ep, 0.0_ep); VR = (0.0_ep, 0.0_ep); select = .true.
        allocate(work(2 * n))
        call ztgevc('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                    mm, mout, work, rwork, info)
        deallocate(work)
        allocate(s_r(n), dif_r(n), s_g(n), dif_g(n))
        allocate(work(max(1, 2 * n * n)), iwork(n+2))
        call ztgsna('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                    s_r, dif_r, mm, mout, work, size(work), iwork, info)
        deallocate(work, iwork)
        call target_ztgsna('B', 'A', select, n, A_s, n, B_s, n, VL, n, VR, n, &
                           s_g, dif_g, mm, mout, info)
        call sort_desc_local(s_r, mout); call sort_desc_local(s_g, mout)
        err = max_rel_err_vec(s_g(1:mout), s_r(1:mout))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, A_s, B_s, Q_s, Z_s, alpha, beta, rwork, &
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
end program test_ztgsna
