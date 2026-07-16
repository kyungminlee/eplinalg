! zgees: complex Schur decomposition. SORT='N'.
program test_zgees
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgees
    use ref_quad_lapack, only: zgees
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, sdim_r, sdim_g, j
    complex(ep), allocatable :: A(:,:), Aref(:,:), Agot(:,:)
    complex(ep), allocatable :: W_r(:), W_g(:), VS_r(:,:), VS_g(:,:), work(:)
    real(ep),    allocatable :: rwork(:), mods_r(:), mods_g(:)
    logical,     allocatable :: bwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zgees', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_complex(n, n, A, seed = 101001 + 47 * i)
        allocate(Aref(n,n), Agot(n,n), W_r(n), W_g(n), VS_r(n,n), VS_g(n,n))
        allocate(rwork(n), bwork(n), mods_r(n), mods_g(n))
        Aref = A; Agot = A
        call zgees('N', 'N', sel_all_c, n, Aref, n, sdim_r, W_r, &
                   VS_r, n, wopt, -1, rwork, bwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zgees('N', 'N', sel_all_c, n, Aref, n, sdim_r, W_r, &
                   VS_r, n, work, lwork, rwork, bwork, info)
        deallocate(work)
        call target_zgees('N', 'N', sel_all_c, n, Agot, n, sdim_g, W_g, &
                          VS_g, n, info)
        do j = 1, n
            mods_r(j) = abs(W_r(j))
            mods_g(j) = abs(W_g(j))
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Aref, Agot, W_r, W_g, VS_r, VS_g, rwork, bwork, mods_r, mods_g)
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
end program test_zgees
