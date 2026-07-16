! Test BLACS xgebs2d / xgebr2d — broadcast of a complex matrix with
! scope='A' (all-grid), 'R' (row), 'C' (column).
program test_xgebs2d
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use test_data,         only: gen_matrix_complex
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col
    use target_blacs,      only: target_name, target_xgebs2d, target_xgebr2d
    use mpi
    implicit none

    integer, parameter :: m = 3, n = 3, lda = 3
    complex(ep), allocatable :: A(:,:), Aref(:,:)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_xgebs2d', target_name, my_rank)

    call gen_matrix_complex(lda, n, Aref, seed=2727)
    allocate(A(lda, n))

    ! ── Scope='A' (all grid) ────────────────────────────────────────
    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_xgebs2d(my_context, 'A', ' ', m, n, A, lda)
        err = 0.0_ep
    else
        A = (0.0_ep, 0.0_ep)
        call target_xgebr2d(my_context, 'A', ' ', m, n, A, lda, 0, 0)
        err = max_rel_err_mat_z(A(1:m, 1:n), Aref(1:m, 1:n))
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
        call target_xgebs2d(my_context, 'R', ' ', m, n, A, lda)
    else
        A = (0.0_ep, 0.0_ep)
        call target_xgebr2d(my_context, 'R', ' ', m, n, A, lda, my_row, 0)
        err = max_rel_err_mat_z(A(1:m, 1:n), Aref(1:m, 1:n))
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
        call target_xgebs2d(my_context, 'C', ' ', m, n, A, lda)
    else
        A = (0.0_ep, 0.0_ep)
        call target_xgebr2d(my_context, 'C', ' ', m, n, A, lda, 0, my_col)
        err = max_rel_err_mat_z(A(1:m, 1:n), Aref(1:m, 1:n))
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
end program test_xgebs2d
