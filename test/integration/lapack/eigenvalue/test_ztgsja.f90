! ztgsja: complex pre-processed GSVD by Jacobi rotations. Smoke test.
! Setup: m=p=n, k=0, l=n with A and B both upper triangular.
program test_ztgsja
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztgsja
    use ref_quad_lapack, only: ztgsja
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, m, p, k, l, info, j, ncyc_r, ncyc_g
    complex(ep), allocatable :: A0(:,:), B0(:,:)
    complex(ep), allocatable :: A_r(:,:), B_r(:,:), U_r(:,:), V_r(:,:), Q_r(:,:)
    complex(ep), allocatable :: A_g(:,:), B_g(:,:), U_g(:,:), V_g(:,:), Q_g(:,:)
    complex(ep), allocatable :: work(:)
    real(ep), allocatable :: alpha_r(:), beta_r(:), alpha_g(:), beta_g(:)
    real(ep) :: tola, tolb
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztgsja', target_name)
    do i = 1, size(ns)
        n = ns(i); m = n; p = n; k = 0; l = n
        call gen_matrix_complex(m, n, A0, seed = 147101 + 47 * i)
        call gen_matrix_complex(p, n, B0, seed = 147111 + 47 * i)
        do j = 1, min(m, n)
            A0(j+1:m, j) = (0.0_ep, 0.0_ep)
            A0(j, j) = cmplx(abs(A0(j, j)) + 1.0_ep, 0.0_ep, ep)
        end do
        do j = 1, min(p, n)
            B0(j+1:p, j) = (0.0_ep, 0.0_ep)
            B0(j, j) = cmplx(abs(B0(j, j)) + 1.0_ep, 0.0_ep, ep)
        end do
        allocate(A_r(m,n), B_r(p,n), U_r(m,m), V_r(p,p), Q_r(n,n))
        allocate(A_g(m,n), B_g(p,n), U_g(m,m), V_g(p,p), Q_g(n,n))
        allocate(alpha_r(n), beta_r(n), alpha_g(n), beta_g(n))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        U_r = (0.0_ep, 0.0_ep); V_r = (0.0_ep, 0.0_ep); Q_r = (0.0_ep, 0.0_ep)
        do j = 1, m; U_r(j,j) = (1.0_ep, 0.0_ep); end do
        do j = 1, p; V_r(j,j) = (1.0_ep, 0.0_ep); end do
        do j = 1, n; Q_r(j,j) = (1.0_ep, 0.0_ep); end do
        U_g = U_r; V_g = V_r; Q_g = Q_r
        tola = real(m * n, ep) * target_eps
        tolb = real(p * n, ep) * target_eps
        allocate(work(2 * n))
        call ztgsja('U', 'V', 'Q', m, p, n, k, l, A_r, m, B_r, p, &
                    tola, tolb, alpha_r, beta_r, U_r, m, V_r, p, Q_r, n, &
                    work, ncyc_r, info)
        deallocate(work)
        call target_ztgsja('U', 'V', 'Q', m, p, n, k, l, A_g, m, B_g, p, &
                           tola, tolb, alpha_g, beta_g, U_g, m, V_g, p, Q_g, n, &
                           ncyc_g, info)
        block
            real(ep) :: r_r(l), r_g(l)
            integer :: kk
            do kk = 1, l
                r_r(kk) = merge(alpha_r(k+kk) / beta_r(k+kk), huge(1.0_ep), &
                                abs(beta_r(k+kk)) > tiny(1.0_ep))
                r_g(kk) = merge(alpha_g(k+kk) / beta_g(k+kk), huge(1.0_ep), &
                                abs(beta_g(k+kk)) > tiny(1.0_ep))
            end do
            call sort_desc_local(r_r, l); call sort_desc_local(r_g, l)
            err = max_rel_err_vec(r_g, r_r)
        end block
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, U_r, V_r, Q_r, A_g, B_g, U_g, V_g, Q_g, &
                   alpha_r, beta_r, alpha_g, beta_g)
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
end program test_ztgsja
