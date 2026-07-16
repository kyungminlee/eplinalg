! Test BLACS qgesd2d / qgerv2d — point-to-point send/recv of a real
! matrix. On a 2x2 grid, rank (0,0) sends an MxN matrix to rank (1,1);
! the receiver verifies byte-equivalent payload (tol = 0).
program test_qgesd2d
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use test_data,         only: gen_matrix_quad
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_qgesd2d, target_qgerv2d
    use mpi
    implicit none

    integer, parameter :: m = 5, n = 4, lda = 7
    real(ep), allocatable :: A(:,:), Aref(:,:)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qgesd2d', target_name, my_rank)

    ! Skip the p2p case when grid is degenerate (1x1) — there's no
    ! distinct receiver. We still report a passing case so the JSON
    ! shape stays uniform across all-mpiexec configurations.
    if (my_nprow < 2 .or. my_npcol < 2) then
        if (my_rank == 0) call report_case('skipped_grid_1x1', 0.0_ep, tol)
        call report_finalize()
        call grid_exit()
        stop
    end if

    ! Reference matrix every rank can compute (same seed).
    call gen_matrix_quad(lda, n, Aref, seed=4242)
    allocate(A(lda, n))
    A = 0.0_ep
    err = 0.0_ep

    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_qgesd2d(my_context, m, n, A, lda, 1, 1)
    else if (my_row == 1 .and. my_col == 1) then
        call target_qgerv2d(my_context, m, n, A, lda, 0, 0)
        ! Compare only the m-by-n payload (lda may exceed m).
        err = max_rel_err_mat(A(1:m, 1:n), Aref(1:m, 1:n))
        ! qgesd2d copies the column-stride (lda) layout exactly, so
        ! the rows m+1..lda inside cols 1..n on the receiver come from
        ! whatever bits MPI moved — those are part of the payload too.
        ! The point of the test is the m×n submatrix; off-payload bits
        ! are not part of the documented contract.
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
end program test_qgesd2d
