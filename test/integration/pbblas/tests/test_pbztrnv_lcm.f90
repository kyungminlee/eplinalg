! Test PBZTRNV — non-square grid exercising LCMP > 1.
!
! Closes the LCMP > 1 half of the "larger / non-square grids" gap in
! tests/pbblas/TODO.md (the LCMQ > 1 half is in test_pbdtrnv_lcm).
!
! Uses grid_init_shape(nprow=1, npcol=4) so LCM(1,4)=4, LCMP=4, LCMQ=1.
! In XDIST='R' the upstream WORK-size formula CEIL(Npb,LCMP)*NB and the
! I=0..LCM-1 dispatch loop only exercise their distinguishing logic
! when LCMP > 1; the existing test_pbztrnv_xdistr runs on a 2×2 grid
! where LCMP=1 and the loop collapses.
!
! Requires connected MPI — see test_pbdtrnv_lcm for the env recipe.
program test_pbztrnv_lcm
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_vec_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_vector_complex
    use pblas_grid,         only: grid_init_shape, grid_exit, my_rank, &
                                  my_context, my_nprow, my_npcol, &
                                  my_row, my_col, numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbztrnv
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: n  = 16
    integer, parameter :: nz = 0
    integer, parameter :: ixrow = 0, ixcol = 0, iyrow = 0, iycol = 0
    integer :: np_loc, nq_loc, lenx, leny, lwork
    integer :: bytes_per_elem, ig, ierr, src_rank, peer
    integer, parameter :: tag_Y = 7641
    complex(ep), allocatable :: X_glob(:), X_loc(:), Y_loc(:)
    complex(ep), allocatable :: Y_glob_got(:), Y_glob_ref(:)
    complex(ep), allocatable :: buf1d(:)
    complex(ep) :: beta
    real(ep)    :: err, tol
    character(len=64) :: label
    integer :: pr, locn, il

    ! 1×4 grid: NPROW=1, NPCOL=4 ⇒ LCMP=4, LCMQ=1.
    call grid_init_shape(1, 4)
    call report_init('pbztrnv_lcm', target_name, my_rank)

    if (my_context < 0) then
        if (my_rank == 0) then
            call report_case('skipped: nproc /= 4 (need connected MPI 4-rank world)', &
                             0.0_ep, 0.0_ep)
        end if
        call report_finalize()
        call grid_exit()
        return
    end if

    ! XDIST='R': X distributed across columns within row IXROW; Y
    ! distributed across rows within column IYCOL.
    nq_loc = numroc_local(n, nb, my_col, ixcol, my_npcol)
    np_loc = numroc_local(n, nb, my_row, iyrow, my_nprow)
    lenx = max(1, nq_loc)
    leny = max(1, np_loc)
    lwork = max(64, 4 * (n + nb))

    call gen_vector_complex(n, X_glob, seed = 9801)
    allocate(X_loc(lenx), Y_loc(leny))
    X_loc = (0.0_ep, 0.0_ep)
    Y_loc = (0.0_ep, 0.0_ep)
    if (my_row == ixrow .and. nq_loc > 0) then
        do ig = 1, n
            block
                integer :: zero_based, blk, owner_c, jl
                zero_based = ig - 1
                blk = zero_based / nb
                owner_c = mod(blk + ixcol, my_npcol)
                if (owner_c == my_col) then
                    jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                    X_loc(jl) = X_glob(ig)
                end if
            end block
        end do
    end if
    beta = (0.0_ep, 0.0_ep)

    call target_pbztrnv(my_context, 'R', 'C', n, nb, nz, X_loc, 1, beta, &
                        Y_loc, 1, ixrow, ixcol, iyrow, iycol, &
                        lenx = lenx, leny = leny, lwork = lwork)

    bytes_per_elem = storage_size((0.0_ep, 0.0_ep)) / 8
    if (my_rank == 0) then
        allocate(Y_glob_got(n))
        Y_glob_got = (0.0_ep, 0.0_ep)
        do pr = 0, my_nprow - 1
            locn = numroc_local(n, nb, pr, iyrow, my_nprow)
            if (locn == 0) cycle
            allocate(buf1d(locn))
            if (pr == my_row .and. my_col == iycol) then
                buf1d = Y_loc(1:locn)
            else
                src_rank = pr * my_npcol + iycol
                call mpi_recv(buf1d, locn * bytes_per_elem, mpi_byte, &
                              src_rank, tag_Y, mpi_comm_world, &
                              mpi_status_ignore, ierr)
            end if
            do ig = 1, n
                block
                    integer :: zero_based, blk, owner_r
                    zero_based = ig - 1
                    blk = zero_based / nb
                    owner_r = mod(blk + iyrow, my_nprow)
                    if (owner_r == pr) then
                        il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                        Y_glob_got(ig) = buf1d(il)
                    end if
                end block
            end do
            deallocate(buf1d)
        end do
    else
        if (my_col == iycol .and. np_loc > 0) then
            peer = 0
            call mpi_send(Y_loc, np_loc * bytes_per_elem, mpi_byte, peer, &
                          tag_Y, mpi_comm_world, ierr)
        end if
    end if

    if (my_rank == 0) then
        allocate(Y_glob_ref(n))
        Y_glob_ref = conjg(X_glob)
        err = max_rel_err_vec_z(Y_glob_got, Y_glob_ref)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0)') 'grid=', my_nprow, 'x', my_npcol, &
                                              ',xdist=R,LCMP=', my_npcol, &
                                              ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(Y_glob_got, Y_glob_ref)
    end if

    deallocate(X_glob, X_loc, Y_loc)

    call report_finalize()
    call grid_exit()
end program test_pbztrnv_lcm
