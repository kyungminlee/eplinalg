! Test BLACS qtrsd2d / qtrrv2d — point-to-point send/recv of a
! trapezoidal real matrix. On a 2x2 grid, rank (0,0) sends the upper
! trapezoid (UPLO='U', DIAG='N') of an MxN matrix to rank (1,1); the
! receiver verifies byte-equivalent payload on the trapezoid (the
! lower triangle is not transmitted; off-trapezoid cells in the
! receive buffer are preserved by the wrapper which preloads At with
! the in-buffer contents before calling the BLACS recv).
!
! BLACS upper-trapezoidal convention (BI_GetMpiTrType.c, M>N branch):
!   col j (1-indexed) carries rows 1..(M-N+j); the extra M-N rows
!   above the j==N diagonal are part of the trapezoid. M==N is the
!   ordinary upper-triangular case.  Test uses M==N to keep the
!   shape unambiguous and exercise the canonical triangular path.
program test_qtrsd2d
    use prec_kinds,        only: ep
    use blacs_prec_report, only: report_init, report_case, report_finalize
    use test_data,         only: gen_matrix_quad
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_row, my_col, my_nprow, my_npcol
    use target_blacs,      only: target_name, target_qtrsd2d, target_qtrrv2d
    use mpi
    implicit none

    integer, parameter :: m = 5, n = 5, lda = 6
    real(ep), allocatable :: A(:,:), Aref(:,:)
    real(ep) :: err
    real(ep), parameter :: tol = 0.0_ep
    real(ep), parameter :: BAD = huge(0.0_ep)
    integer :: i, j, fail_local, fail_global, ierr

    call grid_init()
    call report_init('blacs_qtrsd2d', target_name, my_rank)

    if (my_nprow < 2 .or. my_npcol < 2) then
        if (my_rank == 0) call report_case('skipped_grid_1x1', 0.0_ep, tol)
        call report_finalize(); call grid_exit(); stop
    end if

    call gen_matrix_quad(lda, n, Aref, seed=8181)
    allocate(A(lda, n))
    A = 0.0_ep
    err = 0.0_ep

    if (my_row == 0 .and. my_col == 0) then
        A = Aref
        call target_qtrsd2d(my_context, 'U', 'N', m, n, A, lda, 1, 1)
    else if (my_row == 1 .and. my_col == 1) then
        ! Receiver: the wrapper preloads At with current A contents,
        ! so off-trapezoid cells round-trip through q2t→t2q (exact
        ! at kind16) and stay zero.
        call target_qtrrv2d(my_context, 'U', 'N', m, n, A, lda, 0, 0)
        ! Verify upper triangle (i <= j) is exact.
        do j = 1, n
            do i = 1, min(j, m)
                if (A(i, j) /= Aref(i, j)) err = BAD
            end do
        end do
        ! Verify strictly-lower entries were not overwritten (still 0).
        do j = 1, n
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
        call report_case('upper_trapezoid', err, tol)
    end if

    call report_finalize()
    call grid_exit()
end program test_qtrsd2d
