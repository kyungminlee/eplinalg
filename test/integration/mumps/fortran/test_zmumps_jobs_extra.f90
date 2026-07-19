! Extra JOB-code coverage for ZMUMPS — mirror of test_dmumps_jobs_extra.

program test_zmumps_jobs_extra
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_end, mumps_default_tol
    use mpi
    implicit none

    integer, parameter :: n = 12
    integer            :: ierr, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_solve(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_zmumps_jobs_extra', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 22001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = mumps_default_tol(n)

    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4')
    call run_job(id, 3, 'job=3')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,4,3,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 1, 'job=1')
    call run_job(id, 5, 'job=5')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,1,5,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4 (round 1)')
    call run_job(id, 3, 'job=3 (round 1)')
    call run_job(id, -4, 'job=-4 free factors')
    id%RHS = q2t_c(b)
    call run_job(id, 2, 'job=2 re-factor')
    call run_job(id, 3, 'job=3 re-solve')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,4,3,-4,2,3,-2', err, tol)
    call mumps_end(id);  deallocate(x_solve)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine run_job(id, job, label)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: job
        character(len=*),   intent(in)    :: label
        id%JOB = job
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,a,a,i0,a,i0)') trim(label), &
                ' failed: INFOG(1)=', id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if
    end subroutine run_job

end program test_zmumps_jobs_extra
