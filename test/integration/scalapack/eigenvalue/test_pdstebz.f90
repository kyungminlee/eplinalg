program test_pdstebz
    ! PDSTEBZ — bisection eigenvalue solver for symmetric tridiagonal.
    ! Compare its RANGE='A' output against serial LAPACK dstebz.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dstebz
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol
    use target_scalapack,  only: target_name, target_eps, target_pdstebz
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer :: i, n, m_got, m_ref, nsplit_got, nsplit_ref, info, lwork, liwork
    real(ep), allocatable :: d(:), e(:), W_got(:), W_ref(:)
    real(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: iblock_got(:), isplit_got(:), iwork(:)
    integer,  allocatable :: iblock_ref(:), isplit_ref(:), iwork_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdstebz', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)

        allocate(d(n), e(n))
        block
            integer :: jj, s
            do jj = 1, n
                s = 22401 + 31 * i + jj
                d(jj) = real(mod(s,       1009), ep) / 1009.0_ep + real(2 * jj, ep)
                e(jj) = real(mod(s + 113, 1019), ep) / 1019.0_ep * 0.1_ep
            end do
            e(n) = 0.0_ep
        end block

        allocate(W_got(n), iblock_got(n), isplit_got(n))

        ! Upstream PDSTEBZ requires LWORK >= max(5*N, 7),
        ! LIWORK >= max(4*N, 14, NPROW*NPCOL). Pad generously.
        lwork  = max(5 * n, 7) + 64
        liwork = max(4 * n, my_nprow * my_npcol) + 64
        allocate(work(lwork), iwork(liwork))

        call target_pdstebz(my_context, 'A', 'B', n, 0.0_ep, 0.0_ep, 0, 0, &
                            0.0_ep, d, e, m_got, nsplit_got, W_got, &
                            iblock_got, isplit_got, work, lwork, &
                            iwork, liwork, info)

        if (my_rank == 0) then
            allocate(W_ref(n), iblock_ref(n), isplit_ref(n), &
                     work_ref(4 * n), iwork_ref(3 * n))
            call dstebz('A', 'B', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, &
                        d, e, m_ref, nsplit_ref, W_ref, iblock_ref, &
                        isplit_ref, work_ref, iwork_ref, info)

            if (m_got == m_ref) then
                err = max_rel_err_vec(W_got(1:m_got), W_ref(1:m_ref))
            else
                err = 1.0_ep
            end if
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',m=', m_got
            call report_case(trim(label), err, tol)
            deallocate(W_ref, iblock_ref, isplit_ref, work_ref, iwork_ref)
        end if
        deallocate(d, e, W_got, iblock_got, isplit_got, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdstebz
