! dggev3: blocked generalized eig (no Schur form output). JOBVL/VR='N'.
program test_dggev3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dggev3
    use ref_quad_lapack, only: dggev3
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j
    real(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    real(ep), allocatable :: ar_r(:), ai_r(:), be_r(:), ar_g(:), ai_g(:), be_g(:)
    real(ep), allocatable :: vl(:,:), vr(:,:), work(:), m_r(:), m_g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dggev3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 126001 + 47 * i)
        call gen_matrix_quad(n, n, B, seed = 126011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(ar_r(n), ai_r(n), be_r(n), ar_g(n), ai_g(n), be_g(n))
        allocate(vl(1,1), vr(1,1), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call dggev3('N', 'N', n, Aref, n, Bref, n, ar_r, ai_r, be_r, &
                    vl, 1, vr, 1, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dggev3('N', 'N', n, Aref, n, Bref, n, ar_r, ai_r, be_r, &
                    vl, 1, vr, 1, work, lwork, info)
        deallocate(work)
        call target_dggev3('N', 'N', n, Agot, n, Bgot, n, ar_g, ai_g, be_g, &
                           vl, 1, vr, 1, info)
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
                   vl, vr, m_r, m_g)
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
end program test_dggev3
