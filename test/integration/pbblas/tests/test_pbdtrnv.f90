! Test PBDTRNV — distributed vector transpose between column-vector
! and row-vector distributions. Vector analog of pbdtran:
!
!   XDIST = 'C': X is a column vector of length N held in column IXCOL,
!                block-cyclic over NPROW with block size NB.
!                Y is a row vector of length N held in row IYROW,
!                block-cyclic over NPCOL.
!
! Reference: Y_ref(i) = X_glob(i) for TRANS='T' (the conjugation case
! is exercised in test_pbztrnv).
!
! Scope note: the upstream WORK-size formula
!   IXCOL != -1: CEIL(Nqb,LCMQ)*NB
! depends on grid LCM bookkeeping that we don't compute here; we use a
! generous upper bound (4*(N+NB)) that envelopes every 1×1 / 2×2 /
! near-square configuration this test runs under. NZ=0 keeps the
! offset bookkeeping out of the way (TODO.md notes the NZ>0 paths
! remain uncovered).
program test_pbdtrnv
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
    integer, parameter :: ixrow = 0, ixcol = 0, iyrow = 0, iycol = 0
    integer :: np_loc, nq_loc, lenx, leny, lwork
    integer :: bytes_per_elem, ig, ierr, src_rank, peer
    integer, parameter :: tag_Y = 7501
    real(ep), allocatable :: X_glob(:), X_loc(:), Y_loc(:), Y0_loc(:)
    real(ep), allocatable :: Y_glob_got(:), Y_glob_ref(:)
    real(ep), allocatable :: buf1d(:)
    real(ep) :: beta, err, tol
    character(len=64) :: label
    integer :: pc, locn, jl

    call grid_init()
    call report_init('pbdtrnv', target_name, my_rank)

    np_loc = numroc_local(n, nb, my_row, ixrow, my_nprow)
    nq_loc = numroc_local(n, nb, my_col, iycol, my_npcol)
    lenx = max(1, np_loc)
    leny = max(1, nq_loc)
    lwork = max(64, 4 * (n + nb))

    ! Generate a deterministic global X on every rank, copy local slice
    ! out for the column processes that own X (col == IXCOL).
    call gen_vector_quad(n, X_glob, seed = 9001)
    allocate(X_loc(lenx), Y_loc(leny), Y0_loc(leny))
    X_loc  = 0.0_ep
    Y_loc  = 0.0_ep
    Y0_loc = 0.0_ep
    if (my_col == ixcol .and. np_loc > 0) then
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
    ! Y starts as zero (beta=0 path); pbdtrnv assigns Y := X^T + beta*Y.
    beta = 0.0_ep

    call target_pbdtrnv(my_context, 'C', 'T', n, nb, nz, X_loc, 1, beta, &
                        Y_loc, 1, ixrow, ixcol, iyrow, iycol, &
                        lenx = lenx, leny = leny, lwork = lwork)

    ! Gather Y_loc → Y_glob_got on rank 0. Y is held in row IYROW,
    ! distributed across all column processes.
    bytes_per_elem = storage_size(0.0_ep) / 8
    if (my_rank == 0) then
        allocate(Y_glob_got(n))
        Y_glob_got = 0.0_ep
        do pc = 0, my_npcol - 1
            locn = numroc_local(n, nb, pc, iycol, my_npcol)
            if (locn == 0) cycle
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
        Y_glob_ref = X_glob   ! transpose of a vector is the vector itself
        err = max_rel_err_vec(Y_glob_got, Y_glob_ref)
        tol = 64.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'xdist=C,n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(Y_glob_got, Y_glob_ref)
    end if

    deallocate(X_glob, X_loc, Y_loc, Y0_loc)

    call report_finalize()
    call grid_exit()
end program test_pbdtrnv
