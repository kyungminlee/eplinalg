! Test BLACS xgamx2d / xgamn2d — complex max/min reduction. For
! complex amx, BLACS uses CABS1 magnitude |re| + |im|.
!
! Each rank loads A(1,1) = (rank, rank) → magnitude 2*rank. Max is
! held by rank nproc-1 (magnitude 2*(nproc-1)); min by rank 0 (mag 0).
program test_xgamx2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nproc, my_npcol
    use target_blacs,      only: target_name, target_xgamx2d, target_xgamn2d
    use mpi
    implicit none

    complex(ep) :: A(1, 1), expected_max, expected_min
    integer  :: rA(1), cA(1)
    integer  :: max_rank, min_rank
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_xgamx2d', target_name, my_rank)

    ! CABS1(rank, rank) = 2*rank → strictly monotone in rank, so the
    ! argmax is rank nproc-1 and argmin is rank 0 (unambiguous).
    max_rank = my_nproc - 1
    min_rank = 0
    expected_max = cmplx(real(max_rank, ep), real(max_rank, ep), ep)
    expected_min = cmplx(real(min_rank, ep), real(min_rank, ep), ep)

    ! Max
    A(1, 1) = cmplx(real(my_rank, ep), real(my_rank, ep), ep)
    rA(1) = -77; cA(1) = -77
    call target_xgamx2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, 1, 0, 0)
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
        call report_case('cabs1_max', err, tol)
    end if

    ! Min
    A(1, 1) = cmplx(real(my_rank, ep), real(my_rank, ep), ep)
    rA(1) = -77; cA(1) = -77
    call target_xgamn2d(my_context, 'A', ' ', 1, 1, A, 1, rA, cA, 1, 0, 0)
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
        call report_case('cabs1_min', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_xgamx2d
