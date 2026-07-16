! zhgeqz: complex QZ algorithm.
program test_zhgeqz
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zhgeqz
    use ref_quad_lapack, only: zgeqrf, zunmqr, zgghrd, zhgeqz
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j
    complex(ep), allocatable :: A(:,:), B(:,:), tau(:), work(:)
    complex(ep), allocatable :: H_r(:,:), T_r(:,:), H_g(:,:), T_g(:,:), Q(:,:), Z(:,:)
    complex(ep), allocatable :: a_r(:), b_r(:), a_g(:), b_g(:)
    real(ep),    allocatable :: m_r(:), m_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhgeqz', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 131001 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 131011 + 47 * i)
        allocate(tau(n))
        call zgeqrf(n, n, B, n, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgeqrf(n, n, B, n, tau, work, lwork, info)
        call zunmqr('L', 'C', n, n, n, B, n, tau, A, n, work, lwork, info)
        deallocate(work)
        do j = 1, n; B(j+1:n, j) = (0.0_ep, 0.0_ep); end do
        allocate(H_r(n,n), T_r(n,n), H_g(n,n), T_g(n,n), Q(n,n), Z(n,n))
        allocate(a_r(n), b_r(n), a_g(n), b_g(n), m_r(n), m_g(n))
        H_r = A; T_r = B; H_g = A; T_g = B
        Q = (0.0_ep, 0.0_ep); Z = (0.0_ep, 0.0_ep)
        do j = 1, n
            Q(j,j) = (1.0_ep, 0.0_ep); Z(j,j) = (1.0_ep, 0.0_ep)
        end do
        call zgghrd('N', 'N', n, 1, n, H_r, n, T_r, n, Q, n, Z, n, info)
        H_g = H_r; T_g = T_r
        call zhgeqz('E', 'N', 'N', n, 1, n, H_r, n, T_r, n, a_r, b_r, &
                    Q, n, Z, n, wopt, -1, m_r, info)  ! reuse m_r as rwork stub
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zhgeqz('E', 'N', 'N', n, 1, n, H_r, n, T_r, n, a_r, b_r, &
                    Q, n, Z, n, work, lwork, m_r, info)
        deallocate(work)
        call target_zhgeqz('E', 'N', 'N', n, 1, n, H_g, n, T_g, n, a_g, b_g, &
                           Q, n, Z, n, info)
        do j = 1, n
            m_r(j) = merge(abs(a_r(j)) / abs(b_r(j)), huge(1.0_ep), &
                           abs(b_r(j)) > tiny(1.0_ep))
            m_g(j) = merge(abs(a_g(j)) / abs(b_g(j)), huge(1.0_ep), &
                           abs(b_g(j)) > tiny(1.0_ep))
        end do
        call sort_desc(m_r, n); call sort_desc(m_g, n)
        err = max_rel_err_vec(m_g, m_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, tau, H_r, T_r, H_g, T_g, Q, Z, &
                   a_r, b_r, a_g, b_g, m_r, m_g)
    end do
    call report_finalize()
contains
    subroutine sort_desc(x, m)
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
end program test_zhgeqz
