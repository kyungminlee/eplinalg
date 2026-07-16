! Test BLACS xgesd2d / xgerv2d — point-to-point send/recv of a complex
! matrix.
program test_xgesd2d
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use test_data,         only: gen_matrix_complex
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_xgesd2d, target_xgerv2d
    use mpi
    implicit none

    integer, parameter :: m = 4, n = 3, lda = 6
    complex(ep), allocatable :: A(:,:), Aref(:,:)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_xgesd2d', target_name, my_rank)

    if (my_nprow < 2 .or. my_npcol < 2) then
        if (my_rank == 0) call report_case('skipped_grid_1x1', 0.0_ep, tol)
        call report_finalize(); call grid_exit(); stop
    end if

    call gen_matrix_complex(lda, n, Aref, seed=5151)
    allocate(A(lda, n))
    A = (0.0_ep, 0.0_ep)
    err = 0.0_ep

    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_xgesd2d(my_context, m, n, A, lda, 1, 1)
    else if (my_row == 1 .and. my_col == 1) then
        call target_xgerv2d(my_context, m, n, A, lda, 0, 0)
        err = max_rel_err_mat_z(A(1:m, 1:n), Aref(1:m, 1:n))
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
end program test_xgesd2d
