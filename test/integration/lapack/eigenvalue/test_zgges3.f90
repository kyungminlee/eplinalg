! zgges3: complex blocked generalized Schur. SORT='N'.
program test_zgges3
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgges3
    use ref_quad_lapack, only: zgges3
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim_r, sdim_g, j
    complex(ep), allocatable :: A(:,:), B(:,:), Aref(:,:), Bref(:,:), Agot(:,:), Bgot(:,:)
    complex(ep), allocatable :: a_r(:), b_r(:), a_g(:), b_g(:), vsl(:,:), vsr(:,:), work(:)
    real(ep),    allocatable :: rwork(:), m_r(:), m_g(:)
    logical,     allocatable :: bwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgges3', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 123001 + 47 * i)
        call gen_matrix_complex(n, n, B, seed = 123011 + 47 * i)
        allocate(Aref(n,n), Bref(n,n), Agot(n,n), Bgot(n,n))
        allocate(a_r(n), b_r(n), a_g(n), b_g(n), vsl(n,n), vsr(n,n))
        allocate(rwork(8*n), bwork(n), m_r(n), m_g(n))
        Aref = A; Agot = A; Bref = B; Bgot = B
        call zgges3('N', 'N', 'N', sel_all_2, n, Aref, n, Bref, n, sdim_r, &
                    a_r, b_r, vsl, n, vsr, n, wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgges3('N', 'N', 'N', sel_all_2, n, Aref, n, Bref, n, sdim_r, &
                    a_r, b_r, vsl, n, vsr, n, work, lwork, rwork, bwork, info)
        deallocate(work)
        call target_zgges3('N', 'N', 'N', sel_all_2, n, Agot, n, Bgot, n, sdim_g, &
                           a_g, b_g, vsl, n, vsr, n, info)
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
        deallocate(Aref, Bref, Agot, Bgot, a_r, b_r, a_g, b_g, vsl, vsr, &
                   rwork, bwork, m_r, m_g)
    end do
    call report_finalize()
contains
    logical function sel_all_2(a, b)
        complex(ep), intent(in) :: a, b
        sel_all_2 = .true.
        if (a == a) continue
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
end program test_zgges3
