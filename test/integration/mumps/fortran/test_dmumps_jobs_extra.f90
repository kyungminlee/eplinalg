! Extra JOB-code coverage for DMUMPS — supplements test_dmumps_jobs.
!   JOB =  4  — analyze + factor (= 1+2)
!   JOB =  5  — factor + solve   (= 2+3)
!   JOB = -4  — free factors but keep analysis state (allows re-factor
!               with same pattern but new numerics; new factor must
!               match the original solve)
!
! JOB = -3 / 7 / 8 are conditional on save/restore being compiled in
! and need on-disk paths configured; deferred. JOB = 9 is the
! "suggest IRHS_LOC distribution" entrypoint; only meaningful when
! the user is preparing a distributed-RHS run, deferred.

program test_dmumps_jobs_extra
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, &
                                     dmumps_struc, target_qmumps, &
                                     q2t_r, t2q_r
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_end, mumps_default_tol
    use mpi
    implicit none

    integer, parameter :: n = 12
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_solve(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_dmumps_jobs_extra', target_name)

    call gen_dense_problem(n, A, x_true, b, seed = 11001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = mumps_default_tol(n)

    ! ── JOB sequence -1 → 4 (analyze+factor) → 3 (solve) → -2 ─────
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4 analyze+factor')
    call run_job(id, 3, 'job=3 solve after job=4')
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('seq=-1,4,3,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    ! ── JOB sequence -1 → 1 (analyze) → 5 (factor+solve) → -2 ─────
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 1, 'job=1 analyze')
    call run_job(id, 5, 'job=5 factor+solve')
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('seq=-1,1,5,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    ! ── JOB sequence -1 → 4 → 3 → -4 (free factors, keep analysis) →
    !    2 (re-factor same A) → 3 (re-solve) → -2.
    !    The post-(-4) re-factor must produce the same solution as the
    !    original — JOB=-4 is supposed to free numerical state only,
    !    not the analysis (elimination tree, ordering).
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4 (round 1)')
    call run_job(id, 3, 'job=3 (round 1)')
    call run_job(id, -4, 'job=-4 free factors')
    ! Replace RHS for round 2 since round 1 overwrote it.
    id%RHS = q2t_r(b)
    call run_job(id, 2, 'job=2 re-factor')
    call run_job(id, 3, 'job=3 re-solve')
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('seq=-1,4,3,-4,2,3,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine run_job(id, job, label)
        type(dmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: job
        character(len=*),   intent(in)    :: label
        id%JOB = job
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,a,a,i0,a,i0)') trim(label), &
                ' failed: INFOG(1)=', id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if
    end subroutine run_job

end program test_dmumps_jobs_extra
