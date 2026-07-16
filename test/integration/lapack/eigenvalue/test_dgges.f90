! dgges: generalized Schur of (A,B). SORT='N' so SELCTG isn't called.
program test_dgges
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgges
    use ref_quad_lapack, only: dgges
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim_r, sdim_g, j
    real(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:), ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: vsl(:,:), vsr(:,:), work(:), m_r(:), m_g(:)
    logical,  allocatable :: bwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgges', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 120001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 120011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(vsl(n,n), vsr(n,n), bwork(n), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call dgges('N', 'N', 'N', sel_all_3, n, Aref, n, Bref, n, sdim_r, &
                   ar_r, ai_r, be_r, vsl, n, vsr, n, wopt, -1, bwork, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgges('N', 'N', 'N', sel_all_3, n, Aref, n, Bref, n, sdim_r, &
                   ar_r, ai_r, be_r, vsl, n, vsr, n, work, lwork, bwork, info)
        deallocate(work)
        call target_dgges('N', 'N', 'N', sel_all_3, n, Agot, n, Bgot, n, sdim_g, &
                          ar_g, ai_g, be_g, vsl, n, vsr, n, info)
        do j = 1, n
            if (abs(be_r(j)) > tiny(1.0_ep)) then
                m_r(j) = sqrt(ar_r(j)**2 + ai_r(j)**2) / abs(be_r(j))
            else
                m_r(j) = huge(1.0_ep)
            end if
            if (abs(be_g(j)) > tiny(1.0_ep)) then
                m_g(j) = sqrt(ar_g(j)**2 + ai_g(j)**2) / abs(be_g(j))
            else
                m_g(j) = huge(1.0_ep)
            end if
        end do
        call sort_desc(m_r, n); call sort_desc(m_g, n)
        err = max_rel_err_vec(m_g, m_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Aref, Bref, Agot, Bgot, ar_r, ai_r, be_r, ar_g, ai_g, be_g, &
                   vsl, vsr, bwork, m_r, m_g)
    end do
    call report_finalize()
contains
    logical function sel_all_3(re, im, b)
        real(ep), intent(in) :: re, im, b
        sel_all_3 = .true.
        if (re == re) continue
        if (im == im) continue
        if (b == b) continue
    end function
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
end program test_dgges
