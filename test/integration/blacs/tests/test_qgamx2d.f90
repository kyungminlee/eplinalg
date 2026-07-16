! Test BLACS qgamx2d / qgamn2d — element-wise max/min reduction with
! optional locator (RA, CA) tracking which (prow, pcol) held the
! extremum.
!
! Real qgamx2d/qgamn2d compare by |value| (see dgamx2d_.c). Each rank
! loads A(1,1) = ((-1)**my_rank) * (my_rank + 1) * 0.5_ep so values
! are -0.5, +1.0, -1.5, +2.0 for nproc=4 (|.| = 0.5, 1.0, 1.5, 2.0).
! After qgamx2d → rank (0,0) holds +2.0 (largest |.|, from rank 3).
! After qgamn2d → rank (0,0) holds -0.5 (smallest |.|, from rank 0).
program test_qgamx2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nproc, my_npcol
    use target_blacs,      only: target_name, target_qgamx2d, target_qgamn2d
    use mpi
    implicit none

    real(ep) :: A(1, 1), expected_max, expected_min, err
    integer  :: rA(1), cA(1)
    integer  :: max_rank, min_rank
    integer  :: ldia
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qgamx2d', target_name, my_rank)

    ! Largest |.| comes from the highest rank (=nproc-1). Sign is
    ! ((-1)**(nproc-1)) * (nproc) * 0.5.
    max_rank = my_nproc - 1
    min_rank = 0
    expected_max = real((-1)**max_rank, ep) * real(max_rank + 1, ep) * 0.5_ep
    expected_min = real((-1)**min_rank, ep) * real(min_rank + 1, ep) * 0.5_ep
    ldia = 1

    ! ── Max with locator ────────────────────────────────────────────
    A(1, 1) = real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep
    rA(1) = -77; cA(1) = -77
    call target_qgamx2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, ldia, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected_max) err = BAD
        ! Row-major grid: rank q -> prow=q/npcol, pcol=q%npcol.
        if (rA(1) /= max_rank / my_npcol) err = BAD
        if (cA(1) /= mod(max_rank, my_npcol)) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('max_with_locator', err, tol)
    end if

    ! ── Min with locator ────────────────────────────────────────────
    A(1, 1) = real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep
    rA(1) = -77; cA(1) = -77
    call target_qgamn2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, ldia, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected_min) err = BAD
        if (rA(1) /= min_rank / my_npcol) err = BAD
        if (cA(1) /= mod(min_rank, my_npcol)) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('min_with_locator', err, tol)
    end if

    ! ── Max without locator (LDIA = -1) ─────────────────────────────
    A(1, 1) = real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep
    call target_qgamx2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, -1, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected_max) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('max_no_locator', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_qgamx2d
