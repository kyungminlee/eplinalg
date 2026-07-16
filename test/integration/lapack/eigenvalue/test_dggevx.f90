! dggevx: generalized eig with optional balancing/condition.
! BALANC='N', JOBVL/VR='N', SENSE='N'.
program test_dggevx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggevx
    use ref_quad_lapack, only: dggevx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, ilo_r, ihi_r, ilo_g, ihi_g, j
    real(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:), ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: vl(:,:), vr(:,:), ls_r(:), rs_r(:), ls_g(:), rs_g(:)
    real(ep), allocatable :: rce_r(:), rcv_r(:), rce_g(:), rcv_g(:), work(:), m_r(:), m_g(:)
    integer,  allocatable :: iwork(:)
    logical,  allocatable :: bwork(:)
    real(ep) :: wopt(1), abnrm_r, bbnrm_r, abnrm_g, bbnrm_g, err, tol
    character(len=48) :: label

    call report_init('dggevx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 128001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 128011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(vl(1,1), vr(1,1), ls_r(n), rs_r(n), ls_g(n), rs_g(n))
        allocate(rce_r(n), rcv_r(n), rce_g(n), rcv_g(n))
        allocate(iwork(n+6), bwork(n), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call dggevx('N', 'N', 'N', 'N', n, Aref, n, Bref, n, ar_r, ai_r, be_r, &
                    vl, 1, vr, 1, ilo_r, ihi_r, ls_r, rs_r, abnrm_r, bbnrm_r, &
                    rce_r, rcv_r, wopt, -1, iwork, bwork, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dggevx('N', 'N', 'N', 'N', n, Aref, n, Bref, n, ar_r, ai_r, be_r, &
                    vl, 1, vr, 1, ilo_r, ihi_r, ls_r, rs_r, abnrm_r, bbnrm_r, &
                    rce_r, rcv_r, work, lwork, iwork, bwork, info)
        deallocate(work)
        call target_dggevx('N', 'N', 'N', 'N', n, Agot, n, Bgot, n, ar_g, ai_g, be_g, &
                           vl, 1, vr, 1, ilo_g, ihi_g, ls_g, rs_g, abnrm_g, bbnrm_g, &
                           rce_g, rcv_g, info)
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
                   vl, vr, ls_r, rs_r, ls_g, rs_g, rce_r, rcv_r, rce_g, rcv_g, &
                   iwork, bwork, m_r, m_g)
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
end program test_dggevx
