! Test BLACS igesd2d / igerv2d — integer point-to-point send/recv.
! Integer BLACS routines are NOT migrated per-target; they always link
! against the original i-prefix entries (no precision conversion).
program test_igesd2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_igesd2d, target_igerv2d
    use mpi
    implicit none

    integer, parameter :: m = 4, n = 3, lda = 4
    integer :: A(lda, n), Aref(lda, n)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: i, j, fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_igesd2d', target_name, my_rank)

    if (my_nprow < 2 .or. my_npcol < 2) then
        if (my_rank == 0) call report_case('skipped_grid_1x1', 0.0_ep, tol)
        call report_finalize(); call grid_exit(); stop
    end if

    do j = 1, n
        do i = 1, lda
            Aref(i, j) = 1000 * j + i
        end do
    end do
    A = 0
    err = 0.0_ep

    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_igesd2d(my_context, m, n, A, lda, 1, 1)
    else if (my_row == 1 .and. my_col == 1) then
        call target_igerv2d(my_context, m, n, A, lda, 0, 0)
        do j = 1, n
            do i = 1, m
                if (A(i, j) /= Aref(i, j)) err = BAD
            end do
        end do
    end if

    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('p2p_00_to_11', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_igesd2d
