! Test BLACS igsum2d, igamx2d, igamn2d — integer reductions.
program test_igsum2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nproc, my_npcol
    use target_blacs,      only: target_name, target_igsum2d, &
                                 target_igamx2d, target_igamn2d
    use mpi
    implicit none

    integer :: A(1, 1), expected_sum, expected_max, expected_min
    integer :: rA(1), cA(1)
    integer :: max_rank, min_rank
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_igsum2d', target_name, my_rank)

    ! Each rank contributes (rank + 1) so values are strictly distinct
    ! (1, 2, …, nproc): argmax = nproc-1, argmin = 0 unambiguously.
    expected_sum = my_nproc * (my_nproc + 1) / 2
    expected_max = my_nproc
    expected_min = 1
    max_rank = my_nproc - 1
    min_rank = 0

    ! Sum
    A(1, 1) = my_rank + 1
    call target_igsum2d(my_context, 'A', ' ', 1, 1, A, 1, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected_sum) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('sum_to_00', err, tol)
    end if

    ! Max
    A(1, 1) = my_rank + 1
    rA(1) = -77; cA(1) = -77
    call target_igamx2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, 1, 0, 0)
    err = 0.0_ep
    if (my_row == 0 .and. my_col == 0) then
        if (A(1, 1) /= expected_max) err = BAD
        if (rA(1) /= max_rank / my_npcol) err = BAD
        if (cA(1) /= mod(max_rank, my_npcol)) err = BAD
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('max_to_00', err, tol)
    end if

    ! Min
    A(1, 1) = my_rank + 1
    rA(1) = -77; cA(1) = -77
    call target_igamn2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, 1, 0, 0)
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
        call report_case('min_to_00', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_igsum2d
