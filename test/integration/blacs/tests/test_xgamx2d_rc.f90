! Test BLACS xgamx2d / xgamn2d (complex) under row ('R') and column ('C')
! scopes. Mirrors test_qgamx2d_rc but with complex values; the BLACS
! amx/amn comparison key is Cabs1 = |Re| + |Im|, so we keep all values
! purely real to make the expected-rank computation match the qgamx2d
! pattern exactly.
program test_xgamx2d_rc
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_xgamx2d, target_xgamn2d
    use mpi
    implicit none

    complex(ep) :: A(1, 1), expected_max, expected_min
    real(ep) :: err
    integer  :: rA(1), cA(1)
    integer  :: rk_max, rk_min, j, q
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_xgamx2d_rc', target_name, my_rank)

    ! ── Row scope: max of cabs1 across this process row ─────────────
    rk_max = my_row * my_npcol
    do j = 1, my_npcol - 1
        q = my_row * my_npcol + j
        if (real(q + 1, ep) > real(rk_max + 1, ep)) rk_max = q
    end do
    expected_max = cmplx(real((-1)**rk_max, ep) * real(rk_max + 1, ep) * 0.5_ep, &
                         0.0_ep, ep)
    A(1, 1) = cmplx(real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep, &
                    0.0_ep, ep)
    rA(1) = -77; cA(1) = -77
    call target_xgamx2d(my_context, 'R', ' ', 1, 1, A, 1, rA, cA, 1, my_row, 0)
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

    ! ── Column scope: min of cabs1 across this process column ───────
    rk_min = my_col
    expected_min = cmplx(real((-1)**rk_min, ep) * real(rk_min + 1, ep) * 0.5_ep, &
                         0.0_ep, ep)
    A(1, 1) = cmplx(real((-1)**my_rank, ep) * real(my_rank + 1, ep) * 0.5_ep, &
                    0.0_ep, ep)
    rA(1) = -77; cA(1) = -77
    call target_xgamn2d(my_context, 'C', ' ', 1, 1, A, 1, rA, cA, 1, 0, my_col)
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
end program test_xgamx2d_rc
