! Test BLACS qgamx2d / qgamn2d under row ('R') and column ('C') scope.
program test_qgamx2d_rc
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_qgamx2d, target_qgamn2d
    use mpi
    implicit none

    real(ep) :: A(1, 1), expected_max, expected_min, err
    integer  :: rA(1), cA(1)
    integer  :: rk_max, rk_min, j, q
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qgamx2d_rc', target_name, my_rank)

    ! Each rank loads ((-1)**rank) * (rank + 1) * 0.5_ep — same data
    ! pattern as test_qgamx2d, just reduced by row / column slice.

    ! ── Row scope: max of |.| across this process row ────────────────
    rk_max = my_row * my_npcol
    do j = 1, my_npcol - 1
        q = my_row * my_npcol + j
        if (real(q + 1, ep) > real(rk_max + 1, ep)) rk_max = q
    end do
    expected_max = real((-1)**rk_max, ep) * real(rk_max + 1, ep) * 0.5_ep
    A(1, 1) = real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep
    rA(1) = -77; cA(1) = -77
    call target_qgamx2d(my_context, 'R', ' ', 1, 1, A, 1, rA, cA, 1, my_row, 0)
    err = 0.0_ep
    if (my_col == 0) then
        if (A(1, 1) /= expected_max) err = BAD
    end if
    fail_local = 0; if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('row_scope_max', err, tol)
    end if

    ! ── Column scope: min of |.| across this process column ─────────
    rk_min = my_col
    expected_min = real((-1)**rk_min, ep) * real(rk_min + 1, ep) * 0.5_ep
    A(1, 1) = real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep
    rA(1) = -77; cA(1) = -77
    call target_qgamn2d(my_context, 'C', ' ', 1, 1, A, 1, rA, cA, 1, 0, my_col)
    err = 0.0_ep
    if (my_row == 0) then
        if (A(1, 1) /= expected_min) err = BAD
    end if
    fail_local = 0; if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('col_scope_min', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_qgamx2d_rc
