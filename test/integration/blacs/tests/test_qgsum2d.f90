! Test BLACS qgsum2d — combine sum reduction. Each rank starts with a
! 1x1 buffer holding its rank-as-real; after qgsum2d on rank (0,0),
! that rank holds 0+1+...+nproc-1 = nproc*(nproc-1)/2.
!
! Also test all-reduce variant (rdest=-1): every rank holds the sum.
program test_qgsum2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nproc
    use target_blacs,      only: target_name, target_qgsum2d
    use mpi
    implicit none

    real(ep) :: A(1, 1), expected, err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qgsum2d', target_name, my_rank)

    expected = real(my_nproc, ep) * real(my_nproc - 1, ep) / 2.0_ep

    ! ── Reduce-to-(0,0) ─────────────────────────────────────────────
    A(1, 1) = real(my_rank, ep)
    call target_qgsum2d(my_context, 'A', ' ', 1, 1, A, 1, 0, 0)
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
    A(1, 1) = real(my_rank, ep)
    call target_qgsum2d(my_context, 'A', ' ', 1, 1, A, 1, -1, -1)
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
end program test_qgsum2d
