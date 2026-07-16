! ztrexc: exchange diagonal elements of a complex Schur form.
program test_ztrexc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_ztrexc
    use ref_quad_lapack, only: zgees, ztrexc
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim, j
    complex(ep), allocatable :: A(:,:), T_r(:,:), T_g(:,:), Q_r(:,:), Q_g(:,:)
    complex(ep), allocatable :: W(:), VS(:,:), work(:)
    real(ep),    allocatable :: rwork(:)
    logical,     allocatable :: bwork(:)
    real(ep),    allocatable :: mods_r(:), mods_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('ztrexc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 107001 + 47 * i)
        allocate(T_r(n,n), T_g(n,n), Q_r(n,n), Q_g(n,n), W(n), VS(n,n))
        allocate(rwork(n), bwork(n))
        T_r = A
        call zgees('V', 'N', sel_all_c, n, T_r, n, sdim, W, VS, n, &
                   wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgees('V', 'N', sel_all_c, n, T_r, n, sdim, W, VS, n, &
                   work, lwork, rwork, bwork, info)
        deallocate(work)
        T_g = T_r; Q_r = VS; Q_g = VS
        call ztrexc('V', n, T_r, n, Q_r, n, n, 1, info)
        call target_ztrexc('V', n, T_g, n, Q_g, n, n, 1, info)
        ! Compare sorted moduli of diagonal entries (eigenvalues).
        allocate(mods_r(n), mods_g(n))
        do j = 1, n
            mods_r(j) = abs(T_r(j, j))
            mods_g(j) = abs(T_g(j, j))
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(T_r, T_g, Q_r, Q_g, W, VS, rwork, bwork, mods_r, mods_g)
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
end program test_ztrexc
