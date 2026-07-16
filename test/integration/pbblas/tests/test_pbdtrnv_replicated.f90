! Test PBDTRNV — replicated input/output paths (IXCOL=-1, IYROW=-1).
!
! Closes the IXCOL=-1 / IYROW=-1 coverage gap previously listed in
! tests/pbblas/TODO.md. Mirror of test_pbdtrnv (XDIST='C', TRANS='T')
! with two changes:
!
!   IXCOL = -1: X is column-distributed and *replicated* across every
!               process column (every column-process holds its own
!               column-cyclic slice of the global X).
!   IYROW = -1: Y is row-distributed and replicated across every
!               process row on output (every row-process holds its
!               own row-cyclic slice of the result).
!
! Reference: Y_ref(i) = X_glob(i) for TRANS='T'.
!
! Verification: gather Y from a single process row (row 0) for the
! "got" result, then independently gather Y from row my_nprow-1 for
! a sanity cross-check that the replicated copies agree byte-for-byte.
program test_pbdtrnv_replicated
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_vec
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_vector_quad
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context, &
                                  my_nprow, my_npcol, my_row, my_col, &
                                  numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbdtrnv
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: n  = 16
    integer, parameter :: nz = 0
    integer, parameter :: ixrow = 0, ixcol = -1, iyrow = -1, iycol = 0
    integer :: np_loc, nq_loc, lenx, leny, lwork
    integer :: bytes_per_elem, ig, ierr, src_rank, peer
    integer, parameter :: tag_Y0 = 7521, tag_YN = 7522
    real(ep), allocatable :: X_glob(:), X_loc(:), Y_loc(:)
    real(ep), allocatable :: Y_glob_got(:), Y_glob_lastrow(:), Y_glob_ref(:)
    real(ep), allocatable :: buf1d(:)
    real(ep) :: beta, err, err_rep, tol
    character(len=64) :: label
    integer :: pc, locn, jl, last_row

    call grid_init()
    call report_init('pbdtrnv_replicated', target_name, my_rank)

    np_loc = numroc_local(n, nb, my_row, ixrow, my_nprow)
    nq_loc = numroc_local(n, nb, my_col, iycol, my_npcol)
    lenx = max(1, np_loc)
    leny = max(1, nq_loc)
    ! Replicated paths use the larger formula:
    !   Size(WORK) = CEIL(Nqb,LCMQ)*NB * MIN(LCMQ,CEIL(NN,NB))
    ! 4*(N+NB) envelopes the 2x2 / 1x1 grid sizes used here.
    lwork = max(64, 4 * (n + nb))

    call gen_vector_quad(n, X_glob, seed = 9301)
    allocate(X_loc(lenx), Y_loc(leny))
    X_loc = 0.0_ep
    Y_loc = 0.0_ep
    ! IXCOL=-1: every process column holds its column-cyclic slice of X.
    if (np_loc > 0) then
        do ig = 1, n
            block
                integer :: zero_based, blk, owner_r, il
                zero_based = ig - 1
                blk = zero_based / nb
                owner_r = mod(blk + ixrow, my_nprow)
                if (owner_r == my_row) then
                    il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                    X_loc(il) = X_glob(ig)
                end if
            end block
        end do
    end if
    beta = 0.0_ep

    call target_pbdtrnv(my_context, 'C', 'T', n, nb, nz, X_loc, 1, beta, &
                        Y_loc, 1, ixrow, ixcol, iyrow, iycol, &
                        lenx = lenx, leny = leny, lwork = lwork)

    ! Gather Y from process row 0 → Y_glob_got, and from row my_nprow-1
    ! → Y_glob_lastrow. With IYROW=-1, both rows must hold identical data.
    bytes_per_elem = storage_size(0.0_ep) / 8
    last_row = my_nprow - 1

    if (my_rank == 0) then
        allocate(Y_glob_got(n), Y_glob_lastrow(n))
        Y_glob_got      = 0.0_ep
        Y_glob_lastrow  = 0.0_ep
        do pc = 0, my_npcol - 1
            locn = numroc_local(n, nb, pc, iycol, my_npcol)
            if (locn == 0) cycle
            allocate(buf1d(locn))

            ! Pull from row 0.
            if (pc == my_col .and. my_row == 0) then
                buf1d = Y_loc(1:locn)
            else
                src_rank = 0 * my_npcol + pc
                call mpi_recv(buf1d, locn * bytes_per_elem, mpi_byte, &
                              src_rank, tag_Y0, mpi_comm_world, &
                              mpi_status_ignore, ierr)
            end if
            do ig = 1, n
                block
                    integer :: zero_based, blk, owner_c
                    zero_based = ig - 1
                    blk = zero_based / nb
                    owner_c = mod(blk + iycol, my_npcol)
                    if (owner_c == pc) then
                        jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                        Y_glob_got(ig) = buf1d(jl)
                    end if
                end block
            end do

            ! Pull from last row (only meaningful when my_nprow > 1).
            if (last_row > 0) then
                if (pc == my_col .and. my_row == last_row) then
                    buf1d = Y_loc(1:locn)
                else
                    src_rank = last_row * my_npcol + pc
                    call mpi_recv(buf1d, locn * bytes_per_elem, mpi_byte, &
                                  src_rank, tag_YN, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    block
                        integer :: zero_based, blk, owner_c
                        zero_based = ig - 1
                        blk = zero_based / nb
                        owner_c = mod(blk + iycol, my_npcol)
                        if (owner_c == pc) then
                            jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                            Y_glob_lastrow(ig) = buf1d(jl)
                        end if
                    end block
                end do
            else
                Y_glob_lastrow = Y_glob_got
            end if

            deallocate(buf1d)
        end do
    else
        if (my_row == 0 .and. nq_loc > 0) then
            peer = 0
            call mpi_send(Y_loc, nq_loc * bytes_per_elem, mpi_byte, peer, &
                          tag_Y0, mpi_comm_world, ierr)
        end if
        if (last_row > 0 .and. my_row == last_row .and. nq_loc > 0) then
            peer = 0
            call mpi_send(Y_loc, nq_loc * bytes_per_elem, mpi_byte, peer, &
                          tag_YN, mpi_comm_world, ierr)
        end if
    end if

    if (my_rank == 0) then
        allocate(Y_glob_ref(n))
        Y_glob_ref = X_glob
        err     = max_rel_err_vec(Y_glob_got, Y_glob_ref)
        err_rep = max_rel_err_vec(Y_glob_lastrow, Y_glob_got)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'replicated(ixcol=-1,iyrow=-1),n=', n, &
                                    ',nb=', nb
        call report_case(trim(label), err, tol)
        ! Replicated copies must match exactly (zero relative error).
        call report_case('replicated_copies_agree', err_rep, 0.0_ep)
        deallocate(Y_glob_got, Y_glob_lastrow, Y_glob_ref)
    end if

    deallocate(X_glob, X_loc, Y_loc)

    call report_finalize()
    call grid_exit()
end program test_pbdtrnv_replicated
