! Test BLACS xgsum2d — complex sum reduction. Each rank starts with a
! 1x1 buffer A(1,1) = (rank, -rank); after qgsum2d the destination
! holds (sum_rank, -sum_rank).
program test_xgsum2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nproc
    use target_blacs,      only: target_name, target_xgsum2d
    use mpi
    implicit none

    complex(ep) :: A(1, 1), expected
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr
    real(ep) :: rsum

    call grid_init()
    call report_init('blacs_xgsum2d', target_name, my_rank)

    rsum = real(my_nproc, ep) * real(my_nproc - 1, ep) / 2.0_ep
    expected = cmplx(rsum, -rsum, ep)

    ! ── Reduce-to-(0,0) ─────────────────────────────────────────────
    A(1, 1) = cmplx(real(my_rank, ep), -real(my_rank, ep), ep)
    call target_xgsum2d(my_context, 'A', ' ', 1, 1, A, 1, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('reduce_to_00', err, tol)
    end if

    ! ── All-reduce (rdest = -1) ─────────────────────────────────────
    A(1, 1) = cmplx(real(my_rank, ep), -real(my_rank, ep), ep)
    call target_xgsum2d(my_context, 'A', ' ', 1, 1, A, 1, -1, -1)
    err = 0.0_ep
    if (A(1, 1) /= expected) err = BAD
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('all_reduce', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_xgsum2d
