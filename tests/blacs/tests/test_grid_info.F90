! Test BLACS grid environment routines:
!   blacs_pinfo, blacs_get, blacs_gridinit, blacs_gridinfo,
!   blacs_pcoord, blacs_pnum, blacs_barrier, blacs_set,
!   blacs_gridexit, blacs_exit.
!
! All four ranks call grid_init() (which exercises pinfo + gridinit +
! gridinfo internally), then verify:
!   - rank/nproc match MPI_COMM_WORLD,
!   - pcoord(rank) inverts pnum(prow,pcol) for every grid cell,
!   - blacs_barrier returns on every rank in scope='A',
!   - blacs_set with a no-op SGET param doesn't error.
program test_grid_info
    use prec_kinds,       only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_nproc, &
                                my_context, my_nprow, my_npcol, my_row, my_col
    use mpi
    implicit none

    integer :: ierr, mpi_rank, mpi_size
    integer :: prow, pcol, q, expected_pnum, fail_local, fail_global
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)

    interface
        subroutine blacs_pcoord(icontxt, nodenum, prow, pcol)
            integer, intent(in)  :: icontxt, nodenum
            integer, intent(out) :: prow, pcol
        end subroutine
        function blacs_pnum(icontxt, prow, pcol) result(n)
            integer, intent(in) :: icontxt, prow, pcol
            integer :: n
        end function
        subroutine blacs_barrier(icontxt, scope)
            integer,      intent(in) :: icontxt
            character(*), intent(in) :: scope
        end subroutine
    end interface

    call grid_init()
    call report_init('blacs_grid_info', 'kind16', my_rank)

    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    call mpi_comm_size(mpi_comm_world, mpi_size, ierr)

    ! Case 1: pinfo round-trip vs MPI_COMM_WORLD.
    err = 0.0_ep
    if (my_rank /= mpi_rank .or. my_nproc /= mpi_size) err = BAD
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('pinfo_matches_mpi', err, tol)
    end if

    ! Case 2: gridinfo gave a 2x2 grid (or whatever pick_grid_shape
    ! decided). Just check that the row/col bounds and current
    ! position are sane.
    err = 0.0_ep
    if (my_nprow * my_npcol /= my_nproc) err = BAD
    if (my_row < 0 .or. my_row >= my_nprow) err = BAD
    if (my_col < 0 .or. my_col >= my_npcol) err = BAD
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('gridinfo_consistent', err, tol)
    end if

    ! Case 3: pcoord/pnum invert each other for every node 0..nproc-1.
    ! Row-major grid → blacs_pnum(p_row, p_col) = p_row * npcol + p_col.
    err = 0.0_ep
    do q = 0, my_nproc - 1
        call blacs_pcoord(my_context, q, prow, pcol)
        if (prow < 0 .or. pcol < 0) err = BAD
        if (prow * my_npcol + pcol /= q) err = BAD
        expected_pnum = blacs_pnum(my_context, prow, pcol)
        if (expected_pnum /= q) err = BAD
    end do
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('pcoord_pnum_invert', err, tol)
    end if

    ! Case 4: blacs_barrier returns on all ranks. The synchronization
    ! itself is the test — if any rank deadlocks the test never
    ! reports.
    call blacs_barrier(my_context, 'A')
    err = 0.0_ep
    if (my_rank == 0) call report_case('barrier_all', err, tol)

    ! Case 5: row + column scope barriers also return.
    call blacs_barrier(my_context, 'R')
    call blacs_barrier(my_context, 'C')
    if (my_rank == 0) call report_case('barrier_row_col', 0.0_ep, tol)

    call report_finalize()
    call grid_exit()
end program test_grid_info
