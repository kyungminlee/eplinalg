! Basic smoke + tolerance baseline for the migrated DMUMPS via the
! kind16 ${LIB_PREFIX}mumps archive (called as `qmumps` here through
! the target_mumps wrapper).
!
! For each n ∈ {8, 32}:
!   1. Generate a quad-precision dense problem A·x_true = b with a
!      diagonally-dominant unsymmetric A (SYM=0).
!   2. Convert A to MUMPS triplet (IRN/JCN/A_trip).
!   3. Solve via the migrated qmumps with JOB=-1 → JOB=6 → JOB=-2.
!   4. Compare the recovered solution to x_true at REAL(KIND=ep).
!
! The test emits per-case JSON via prec_report. Tolerance is the
! repo-standard n³ scaling (matches tests/lapack/linear_solve/test_dgesv.f90).

program test_dmumps_basic
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, dmumps_struc
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_solve6, mumps_check_solution, &
                                     mumps_end
    use mpi
    implicit none

    ! n=1 covers the degenerate single-node elimination tree (no
    ! merges, no panel updates). n=8 / n=32 give realistic small /
    ! medium baselines.
    integer, parameter :: ns(*) = [1, 8, 32]
    integer            :: ierr, i, n, nz, myid
    logical            :: is_host
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc)    :: id
    character(len=48)     :: label

    call MPI_INIT(ierr)
    ! Centralized-input MUMPS: the assembled matrix (N/NNZ/IRN/JCN/A) and
    ! the RHS are supplied on the host only, and the centralized solution
    ! comes back in id%RHS on the host only — every other rank's RHS is
    ! left untouched. So all matrix setup, the solution check and the JSON
    ! report are host-guarded (mirrors MUMPS's own dsimpletest.F, which
    ! wraps its matrix I/O in ``IF (MYID == 0)``). The MUMPS calls
    ! themselves stay collective on every rank. At np=1 the host is the
    ! only rank, so this is behaviourally identical to the single-rank
    ! path CI exercises.
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    is_host = (myid == 0)
    if (is_host) call report_init('test_dmumps_basic', target_name)

    do i = 1, size(ns)
        n = ns(i)
        if (is_host) then
            call gen_dense_problem(n, A, x_true, b, seed = 1001 + 31 * i)
            call dense_to_triplet(A, irn, jcn, A_trip, nz)
        end if

        ! Initialize MUMPS instance (collective).
        call mumps_begin(id, MPI_COMM_WORLD, 0)

        ! Populate problem data on the host only (centralized input).
        if (is_host) call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)

        ! Combined analysis + factorization + solve (collective).
        call mumps_solve6(id)

        ! On exit, id%RHS holds the centralized solution on the host.
        if (is_host) then
            write(label, '(a,i0)') 'n=', n
            call mumps_check_solution(id, x_true, trim(label))
        end if

        ! Destroy MUMPS instance (collective).
        call mumps_end(id)

        if (is_host) deallocate(A, x_true, b, irn, jcn, A_trip)
    end do

    if (is_host) call report_finalize()
    call MPI_FINALIZE(ierr)
    if (is_host) call report_check_status()
end program test_dmumps_basic
