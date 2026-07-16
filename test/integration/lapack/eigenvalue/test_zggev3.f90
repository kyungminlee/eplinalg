! zggev3: complex blocked generalized eig. JOBVL/VR='N'.
program test_zggev3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zggev3
    use ref_quad_lapack, only: zggev3
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, j
    complex(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    complex(ep), allocatable :: a_r(:), b_r(:), a_g(:), b_g(:), vl(:,:), vr(:,:), work(:)
    real(ep),    allocatable :: rwork(:), m_r(:), m_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zggev3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 127001 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 127011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(a_r(n), b_r(n), a_g(n), b_g(n), vl(1,1), vr(1,1))
        allocate(rwork(8*n), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call zggev3('N', 'N', n, Aref, n, Bref, n, a_r, b_r, &
                    vl, 1, vr, 1, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zggev3('N', 'N', n, Aref, n, Bref, n, a_r, b_r, &
                    vl, 1, vr, 1, work, lwork, rwork, info)
        deallocate(work)
        call target_zggev3('N', 'N', n, Agot, n, Bgot, n, a_g, b_g, &
                           vl, 1, vr, 1, info)
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
        deallocate(Aref, Bref, Agot, Bgot, a_r, b_r, a_g, b_g, vl, vr, &
                   rwork, m_r, m_g)
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
end program test_zggev3
