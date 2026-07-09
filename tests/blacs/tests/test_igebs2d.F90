! Test BLACS igebs2d / igebr2d — integer broadcast with scope='A',
! 'R', 'C'.
program test_igebs2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col
    use target_blacs,      only: target_name, target_igebs2d, target_igebr2d
    use mpi
    implicit none

    integer, parameter :: m = 3, n = 2, lda = 3
    integer :: A(lda, n), Aref(lda, n)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: i, j, fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_igebs2d', target_name, my_rank)

    do j = 1, n
        do i = 1, lda
            Aref(i, j) = 100 * j + i + 7
        end do
    end do

    ! ── Scope='A' (all grid) ────────────────────────────────────────
    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_igebs2d(my_context, 'A', ' ', m, n, A, lda)
        err = 0.0_ep
    else
        A = 0
        call target_igebr2d(my_context, 'A', ' ', m, n, A, lda, 0, 0)
        err = 0.0_ep
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
        call report_case('scope_all', err, tol)
    end if

    ! ── Scope='R' (row) — bcast from (my_row, 0) within each row. ───
    err = 0.0_ep
    if (my_col == 0) then
        A = Aref
        call target_igebs2d(my_context, 'R', ' ', m, n, A, lda)
    else
        A = 0
        call target_igebr2d(my_context, 'R', ' ', m, n, A, lda, my_row, 0)
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
        call report_case('scope_row', err, tol)
    end if

    ! ── Scope='C' (column) — bcast from (0, my_col) within each col.─
    err = 0.0_ep
    if (my_row == 0) then
        A = Aref
        call target_igebs2d(my_context, 'C', ' ', m, n, A, lda)
    else
        A = 0
        call target_igebr2d(my_context, 'C', ' ', m, n, A, lda, 0, my_col)
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
        call report_case('scope_col', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_igebs2d
