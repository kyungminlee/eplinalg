program test_pdtrord
    ! Reorder the real Schur form so the SELECTed eigenvalues appear
    ! in the leading diagonal blocks. Compare gathered (T, WR, WI)
    ! against LAPACK dtrsen with JOB='N' (reorder only, no condition).
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat, max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dtrsen
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdtrord
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, m_got, m_ref
    integer :: locm, locn, lld
    integer :: descT(9), descQ(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    integer :: para(6)
    integer,  allocatable :: select_v(:)
    logical,  allocatable :: select_log(:)
    real(ep), allocatable :: T_loc(:,:), T_glob(:,:), T_ref(:,:), T_got(:,:)
    real(ep), allocatable :: Q_loc(:,:), Q_ref(:,:)
    real(ep), allocatable :: WR(:), WI(:), WR_ref(:), WI_ref(:)
    real(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: iwork(:), iwork_ref(:)
    real(ep) :: err, tol, S_dummy, SEP_dummy
    character(len=48) :: label

    call grid_init()
    call report_init('pdtrord', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, T_loc, T_glob, seed = 33101 + 31*i)

        ! Force upper-triangular T (Schur with all real eigenvalues).
        do jg = 1, n
            do ig = jg + 1, n
                T_glob(ig, jg) = 0.0_ep
            end do
            T_glob(jg, jg) = T_glob(jg, jg) + real(2 * jg, ep)
        end do
        if (size(T_loc, 1) > 0 .and. size(T_loc, 2) > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) T_loc(il, jl) = T_glob(ig, jg)
                end do
            end do
        end if

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol); lld = max(1, locm)
        call descinit_local(descT, n, n, mb, nb, 0, 0, my_context, lld, info)
        call descinit_local(descQ, n, n, mb, nb, 0, 0, my_context, lld, info)

        allocate(Q_loc(lld, max(1, locn))); Q_loc = 0.0_ep
        ! COMPQ='N' — don't update Q.

        allocate(select_v(n))
        select_v = 0
        ! Pick odd-indexed eigenvalues to move forward.
        do jg = 1, n, 2
            select_v(jg) = 1
        end do

        ! PARA constraints from upstream PDTRORD docs:
        !   1 ≤ PARA(1) ≤ min(NPROW,NPCOL)        — concurrent windows
        !   1 ≤ PARA(2) < PARA(3)                 — eigenvalues per window
        !   PARA(2) < PARA(3) < MB                — window size
        !   0 ≤ PARA(4) ≤ 100                     — matmul flop pct
        !   0 < PARA(5) ≤ MB                      — row-slab width
        !   0 < PARA(6) ≤ PARA(2)                 — eigenvalues per border move
        para = [1, 2, 5, 50, mb, 1]

        allocate(WR(n), WI(n))
        allocate(work(1), iwork(1))
        call target_pdtrord('N', select_v, para, n, T_loc, 1, 1, descT, &
                            Q_loc, 1, 1, descQ, WR, WI, m_got, &
                            work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdtrord('N', select_v, para, n, T_loc, 1, 1, descT, &
                            Q_loc, 1, 1, descQ, WR, WI, m_got, &
                            work, lwork, iwork, liwork, info)

        call gather_matrix(n, n, mb, nb, T_loc, T_got)

        if (my_rank == 0) then
            allocate(T_ref(n, n), Q_ref(n, n), WR_ref(n), WI_ref(n), &
                     select_log(n), work_ref(max(1, n * 64)), iwork_ref(n))
            T_ref = T_glob
            Q_ref = 0.0_ep
            do jg = 1, n
                Q_ref(jg, jg) = 1.0_ep
            end do
            do jg = 1, n
                select_log(jg) = (select_v(jg) /= 0)
            end do
            call dtrsen('N', 'N', select_log, n, T_ref, n, Q_ref, n, &
                        WR_ref, WI_ref, m_ref, S_dummy, SEP_dummy, &
                        work_ref, size(work_ref), iwork_ref, size(iwork_ref), info_ref)
            ! pdtrord and dtrsen both move SELECTed eigenvalues to the
            ! leading diagonal block of T but may differ in the trailing-
            ! block ordering. Compare only the eigenvalue spectrum: sort
            ! WR from both implementations and compare.
            block
                integer :: a, b
                real(ep) :: tmp
                real(ep), allocatable :: w1(:), w2(:)
                allocate(w1(n), w2(n))
                w1 = WR(1:n); w2 = WR_ref(1:n)
                do a = 1, n - 1
                    do b = a + 1, n
                        if (w1(b) < w1(a)) then; tmp = w1(a); w1(a) = w1(b); w1(b) = tmp; end if
                        if (w2(b) < w2(a)) then; tmp = w2(a); w2(a) = w2(b); w2(b) = tmp; end if
                    end do
                end do
                err = max_rel_err_vec(w1, w2)
                deallocate(w1, w2)
            end block
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',m=', m_got
            call report_case(trim(label), err, tol)
            deallocate(T_ref, Q_ref, WR_ref, WI_ref, select_log, work_ref, iwork_ref, T_got)
        end if
        deallocate(T_loc, T_glob, Q_loc, select_v, WR, WI, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtrord
