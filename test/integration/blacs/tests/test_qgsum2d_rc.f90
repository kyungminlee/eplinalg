! Test BLACS qgsum2d under row ('R') and column ('C') scope. Under a
! real 2x2 grid each row sums independently to its row's destination
! and similarly for columns; in the singleton sandbox the scope
! degenerates to local and the rank holds its own value back.
program test_qgsum2d_rc
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_qgsum2d
    use mpi
    implicit none

    real(ep) :: A(1, 1), expected_row, expected_col, err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr, j

    call grid_init()
    call report_init('blacs_qgsum2d_rc', target_name, my_rank)

    ! ── Row scope: sum of all ranks in this process row to col 0 ─────
    expected_row = 0.0_ep
    do j = 0, my_npcol - 1
        expected_row = expected_row + real(my_row * my_npcol + j, ep)
    end do
    A(1, 1) = real(my_rank, ep)
    call target_qgsum2d(my_context, 'R', ' ', 1, 1, A, 1, my_row, 0)
    err = 0.0_ep
    if (my_col == 0) then
        if (A(1, 1) /= expected_row) err = BAD
    end if
    fail_local = 0; if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('row_scope_to_col0', err, tol)
    end if

    ! ── Column scope: sum of all ranks in this process column to row 0
    expected_col = 0.0_ep
    do j = 0, my_nprow - 1
        expected_col = expected_col + real(j * my_npcol + my_col, ep)
    end do
    A(1, 1) = real(my_rank, ep)
    call target_qgsum2d(my_context, 'C', ' ', 1, 1, A, 1, 0, my_col)
    err = 0.0_ep
    if (my_row == 0) then
        if (A(1, 1) /= expected_col) err = BAD
    end if
    fail_local = 0; if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('col_scope_to_row0', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_qgsum2d_rc
