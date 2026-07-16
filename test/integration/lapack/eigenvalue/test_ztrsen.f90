! ztrsen: complex Schur reordering with SELECT. JOB='N', COMPQ='N'.
program test_ztrsen
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrsen
    use ref_quad_lapack, only: zgees, ztrsen
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim, m_r, m_g, j
    complex(ep), allocatable :: A(:,:), T_r(:,:), T_g(:,:), Q(:,:)
    complex(ep), allocatable :: W_r(:), W_g(:), VS(:,:), work(:)
    real(ep),    allocatable :: rwork(:), mods_r(:), mods_g(:)
    logical,     allocatable :: bwork(:), select(:)
    complex(ep) :: wopt(1)
    real(ep) :: s_r, sep_r, s_g, sep_g, err, tol
    character(len=48) :: label

    call report_init('ztrsen', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 109001 + 47 * i)
        allocate(T_r(n,n), T_g(n,n), Q(1,1), W_r(n), W_g(n), VS(n,n))
        allocate(rwork(n), bwork(n), select(n))
        T_r = A
        call zgees('N', 'N', sel_all_c, n, T_r, n, sdim, W_r, VS, n, &
                   wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgees('N', 'N', sel_all_c, n, T_r, n, sdim, W_r, VS, n, &
                   work, lwork, rwork, bwork, info)
        deallocate(work)
        T_g = T_r
        select = .false.
        do j = 1, n / 2
            select(j) = .true.
        end do
        call ztrsen('N', 'N', select, n, T_r, n, Q, 1, W_r, m_r, &
                    s_r, sep_r, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call ztrsen('N', 'N', select, n, T_r, n, Q, 1, W_r, m_r, &
                    s_r, sep_r, work, lwork, info)
        deallocate(work)
        call target_ztrsen('N', 'N', select, n, T_g, n, Q, 1, W_g, m_g, &
                           s_g, sep_g, info)
        allocate(mods_r(n), mods_g(n))
        do j = 1, n
            mods_r(j) = abs(W_r(j))
            mods_g(j) = abs(W_g(j))
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(T_r, T_g, Q, W_r, W_g, VS, rwork, bwork, select, &
                   mods_r, mods_g)
    end do
    call report_finalize()
contains
    logical function sel_all_c(z)
        complex(ep), intent(in) :: z
        sel_all_c = .true.
        if (z == z) continue
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
end program test_ztrsen
