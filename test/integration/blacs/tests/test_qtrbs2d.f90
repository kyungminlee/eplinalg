! Test BLACS qtrbs2d / qtrbr2d — trapezoidal broadcast with scope='A'
! (all-grid), 'R' (row), 'C' (column). For each scope, the broadcaster
! sends the upper trapezoid of a known matrix; receivers in scope
! verify the trapezoid is exact and off-trapezoid cells stay zero.
program test_qtrbs2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use test_data,         only: gen_matrix_quad
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col
    use target_blacs,      only: target_name, target_qtrbs2d, target_qtrbr2d
    use mpi
    implicit none

    integer, parameter :: m = 5, n = 5, lda = 5
    real(ep), allocatable :: A(:,:), Aref(:,:)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: i, j, fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qtrbs2d', target_name, my_rank)

    call gen_matrix_quad(lda, n, Aref, seed=909)
    allocate(A(lda, n))

    ! ── Scope='A' (all grid) ────────────────────────────────────────
    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_qtrbs2d(my_context, 'A', ' ', 'U', 'N', m, n, A, lda)
        err = 0.0_ep
    else
        A = 0.0_ep
        call target_qtrbr2d(my_context, 'A', ' ', 'U', 'N', m, n, A, lda, 0, 0)
        err = 0.0_ep
        do j = 1, n
            do i = 1, min(j, m)
                if (A(i, j) /= Aref(i, j)) err = BAD
            end do
            do i = j + 1, m
                if (A(i, j) /= 0.0_ep) err = BAD
            end do
        end do
    end if
    fail_local = 0
    if (err > tol) fail_local = 1
    call mpi_allreduce(fail_local, fail_global, 1, mpi_integer, &
                       mpi_max, mpi_comm_world, ierr)
    if (my_rank == 0) then
        if (fail_global /= 0 .and. err <= tol) err = BAD
        call report_case('upper_trapezoid_bcast', err, tol)
    end if

    ! ── Scope='R' (row) — bcast from (my_row, 0) within each row. ───
    err = 0.0_ep
    if (my_col == 0) then
        A = Aref
        call target_qtrbs2d(my_context, 'R', ' ', 'U', 'N', m, n, A, lda)
    else
        A = 0.0_ep
        call target_qtrbr2d(my_context, 'R', ' ', 'U', 'N', m, n, A, lda, my_row, 0)
        do j = 1, n
            do i = 1, min(j, m)
                if (A(i, j) /= Aref(i, j)) err = BAD
            end do
            do i = j + 1, m
                if (A(i, j) /= 0.0_ep) err = BAD
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
        call target_qtrbs2d(my_context, 'C', ' ', 'U', 'N', m, n, A, lda)
    else
        A = 0.0_ep
        call target_qtrbr2d(my_context, 'C', ' ', 'U', 'N', m, n, A, lda, 0, my_col)
        do j = 1, n
            do i = 1, min(j, m)
                if (A(i, j) /= Aref(i, j)) err = BAD
            end do
            do i = j + 1, m
                if (A(i, j) /= 0.0_ep) err = BAD
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
end program test_qtrbs2d
