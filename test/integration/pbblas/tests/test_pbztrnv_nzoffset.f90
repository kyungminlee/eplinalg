! Test PBZTRNV — NZ > 0 (block-offset start). Complex mirror of
! test_pbdtrnv_nzoffset with TRANS='C' (Y_ref = conjg(X_glob)). See
! that file for full commentary on the NZ-aware local layout.
program test_pbztrnv_nzoffset
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_vec_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_vector_complex
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context, &
                                  my_nprow, my_npcol, my_row, my_col, &
                                  numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbztrnv
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: n  = 12
    integer, parameter :: nz = 2
    integer, parameter :: ixrow = 0, ixcol = 0, iyrow = 0, iycol = 0
    integer, parameter :: nn = n + nz
    integer :: np_loc, nq_loc, lenx, leny, lwork
    integer :: bytes_per_elem, ig, ierr, src_rank, peer
    integer, parameter :: tag_Y = 7631
    complex(ep), allocatable :: X_glob(:), X_loc(:), Y_loc(:)
    complex(ep), allocatable :: Y_glob_got(:), Y_glob_ref(:)
    complex(ep), allocatable :: buf1d(:)
    complex(ep) :: beta
    real(ep)    :: err, tol
    character(len=64) :: label
    integer :: pc, locn

    call grid_init()
    call report_init('pbztrnv_nzoffset', target_name, my_rank)

    np_loc = numroc_local(nn, nb, my_row, ixrow, my_nprow)
    if (my_row == ixrow) np_loc = np_loc - nz
    nq_loc = numroc_local(nn, nb, my_col, iycol, my_npcol)
    if (my_col == iycol) nq_loc = nq_loc - nz
    lenx = max(1, np_loc)
    leny = max(1, nq_loc)
    lwork = max(64, 4 * (nn + nb))

    call gen_vector_complex(n, X_glob, seed = 9601)
    allocate(X_loc(lenx), Y_loc(leny))
    X_loc = (0.0_ep, 0.0_ep)
    Y_loc = (0.0_ep, 0.0_ep)
    if (my_col == ixcol .and. np_loc > 0) then
        do ig = 1, n
            block
                integer :: ext_idx, blk, owner_r, cyc, within, il
                ext_idx = nz + ig - 1
                blk = ext_idx / nb
                owner_r = mod(blk + ixrow, my_nprow)
                if (owner_r == my_row) then
                    cyc = blk / my_nprow
                    within = mod(ext_idx, nb)
                    if (my_row == ixrow) then
                        if (cyc == 0) then
                            il = within - nz + 1
                        else
                            il = (nb - nz) + (cyc - 1) * nb + within + 1
                        end if
                    else
                        il = cyc * nb + within + 1
                    end if
                    X_loc(il) = X_glob(ig)
                end if
            end block
        end do
    end if
    beta = (0.0_ep, 0.0_ep)

    call target_pbztrnv(my_context, 'C', 'C', n, nb, nz, X_loc, 1, beta, &
                        Y_loc, 1, ixrow, ixcol, iyrow, iycol, &
                        lenx = lenx, leny = leny, lwork = lwork)

    bytes_per_elem = storage_size((0.0_ep, 0.0_ep)) / 8
    if (my_rank == 0) then
        allocate(Y_glob_got(n))
        Y_glob_got = (0.0_ep, 0.0_ep)
        do pc = 0, my_npcol - 1
            locn = numroc_local(nn, nb, pc, iycol, my_npcol)
            if (pc == iycol) locn = locn - nz
            if (locn <= 0) cycle
            allocate(buf1d(locn))
            if (pc == my_col .and. my_row == iyrow) then
                buf1d = Y_loc(1:locn)
            else
                src_rank = iyrow * my_npcol + pc
                call mpi_recv(buf1d, locn * bytes_per_elem, mpi_byte, &
                              src_rank, tag_Y, mpi_comm_world, &
                              mpi_status_ignore, ierr)
            end if
            do ig = 1, n
                block
                    integer :: ext_idx, blk, owner_c, cyc, within, jl
                    ext_idx = nz + ig - 1
                    blk = ext_idx / nb
                    owner_c = mod(blk + iycol, my_npcol)
                    if (owner_c == pc) then
                        cyc = blk / my_npcol
                        within = mod(ext_idx, nb)
                        if (pc == iycol) then
                            if (cyc == 0) then
                                jl = within - nz + 1
                            else
                                jl = (nb - nz) + (cyc - 1) * nb + within + 1
                            end if
                        else
                            jl = cyc * nb + within + 1
                        end if
                        Y_glob_got(ig) = buf1d(jl)
                    end if
                end block
            end do
            deallocate(buf1d)
        end do
    else
        if (my_row == iyrow .and. nq_loc > 0) then
            peer = 0
            call mpi_send(Y_loc, nq_loc * bytes_per_elem, mpi_byte, peer, &
                          tag_Y, mpi_comm_world, ierr)
        end if
    end if

    if (my_rank == 0) then
        allocate(Y_glob_ref(n))
        Y_glob_ref = conjg(X_glob)
        err = max_rel_err_vec_z(Y_glob_got, Y_glob_ref)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'xdist=C,nz=', nz, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(Y_glob_got, Y_glob_ref)
    end if

    deallocate(X_glob, X_loc, Y_loc)

    call report_finalize()
    call grid_exit()
end program test_pbztrnv_nzoffset
