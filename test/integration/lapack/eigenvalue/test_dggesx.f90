! dggesx: generalized Schur with sensitivity. SORT='N', SENSE='N'.
program test_dggesx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggesx
    use ref_quad_lapack, only: dggesx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, liwork, sdim_r, sdim_g, j, iwopt(1)
    real(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:), ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: vsl(:,:), vsr(:,:), work(:), m_r(:), m_g(:)
    integer,  allocatable :: iwork(:)
    logical,  allocatable :: bwork(:)
    real(ep) :: wopt(1), rce_r(2), rcv_r(2), rce_g(2), rcv_g(2), err, tol
    character(len=48) :: label

    call report_init('dggesx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 124001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 124011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(vsl(n,n), vsr(n,n), bwork(n), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call dggesx('N', 'N', 'N', sel_all_3, 'N', n, Aref, n, Bref, n, sdim_r, &
                    ar_r, ai_r, be_r, vsl, n, vsr, n, rce_r, rcv_r, &
                    wopt, -1, iwopt, -1, bwork, info)
        lwork = max(1, int(wopt(1))); liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dggesx('N', 'N', 'N', sel_all_3, 'N', n, Aref, n, Bref, n, sdim_r, &
                    ar_r, ai_r, be_r, vsl, n, vsr, n, rce_r, rcv_r, &
                    work, lwork, iwork, liwork, bwork, info)
        deallocate(work, iwork)
        call target_dggesx('N', 'N', 'N', sel_all_3, 'N', n, Agot, n, Bgot, n, sdim_g, &
                           ar_g, ai_g, be_g, vsl, n, vsr, n, rce_g, rcv_g, info)
        do j = 1, n
            m_r(j) = merge(sqrt(ar_r(j)**2 + ai_r(j)**2) / abs(be_r(j)), &
                           huge(1.0_ep), abs(be_r(j)) > tiny(1.0_ep))
            m_g(j) = merge(sqrt(ar_g(j)**2 + ai_g(j)**2) / abs(be_g(j)), &
                           huge(1.0_ep), abs(be_g(j)) > tiny(1.0_ep))
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
end program test_dggesx
