! dgeesx: extended Schur decomposition with sensitivity. SORT='N',
! SENSE='N' (no condition estimate exercised — keeps the test small).
program test_dgeesx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgeesx
    use ref_quad_lapack, only: dgeesx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, info, lwork, liwork, sdim_r, sdim_g, j, iwopt(1)
    real(ep), allocatable :: A(:,:), Aref(:,:), Agot(:,:)
    real(ep), allocatable :: WR_r(:), WI_r(:), WR_g(:), WI_g(:)
    real(ep), allocatable :: VS_r(:,:), VS_g(:,:), work(:), mods_r(:), mods_g(:)
    integer,  allocatable :: iwork(:)
    logical,  allocatable :: bwork(:)
    real(ep) :: wopt(1), rce_r, rcv_r, rce_g, rcv_g, err, tol
    character(len=48) :: label

    call report_init('dgeesx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_matrix_quad(n, n, A, seed = 102001 + 47 * i)
        allocate(Aref(n,n), Agot(n,n), WR_r(n), WI_r(n), WR_g(n), WI_g(n))
        allocate(VS_r(n,n), VS_g(n,n), bwork(n), mods_r(n), mods_g(n))
        Aref = A; Agot = A
        call dgeesx('N', 'N', sel_all_r, 'N', n, Aref, n, sdim_r, &
                    WR_r, WI_r, VS_r, n, rce_r, rcv_r, &
                    wopt, -1, iwopt, -1, bwork, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dgeesx('N', 'N', sel_all_r, 'N', n, Aref, n, sdim_r, &
                    WR_r, WI_r, VS_r, n, rce_r, rcv_r, &
                    work, lwork, iwork, liwork, bwork, info)
        deallocate(work, iwork)
        call target_dgeesx('N', 'N', sel_all_r, 'N', n, Agot, n, sdim_g, &
                           WR_g, WI_g, VS_g, n, rce_g, rcv_g, info)
        do j = 1, n
            mods_r(j) = sqrt(WR_r(j)**2 + WI_r(j)**2)
            mods_g(j) = sqrt(WR_g(j)**2 + WI_g(j)**2)
        end do
        call sort_desc(mods_r, n); call sort_desc(mods_g, n)
        err = max_rel_err_vec(mods_g, mods_r)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(Aref, Agot, WR_r, WI_r, WR_g, WI_g, VS_r, VS_g, &
                   bwork, mods_r, mods_g)
    end do
    call report_finalize()
contains
    logical function sel_all_r(re, im)
        real(ep), intent(in) :: re, im
        sel_all_r = .true.
        if (re == re) continue
        if (im == im) continue
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
end program test_dgeesx
