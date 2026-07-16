! Extra JOB-code coverage for ZMUMPS — mirror of test_dmumps_jobs_extra.

program test_zmumps_jobs_extra
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
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
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    call init_id(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4')
    call run_job(id, 3, 'job=3')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,4,3,-2', err, tol)
    call end_id(id);  deallocate(x_solve)

    call init_id(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 1, 'job=1')
    call run_job(id, 5, 'job=5')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,1,5,-2', err, tol)
    call end_id(id);  deallocate(x_solve)

    call init_id(id, n, nz, irn, jcn, A_trip, b)
    call run_job(id, 4, 'job=4 (round 1)')
    call run_job(id, 3, 'job=3 (round 1)')
    call run_job(id, -4, 'job=-4 free factors')
    id%RHS = q2t_c(b)
    call run_job(id, 2, 'job=2 re-factor')
    call run_job(id, 3, 'job=3 re-solve')
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('seq=-1,4,3,-4,2,3,-2', err, tol)
    call end_id(id);  deallocate(x_solve)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine init_id(id, n, nz, irn, jcn, A_trip, b)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        complex(ep),        intent(in)    :: A_trip(:), b(:)
        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_xmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));   id%RHS = q2t_c(b)
    end subroutine init_id

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

    subroutine end_id(id)
        type(zmumps_struc), intent(inout) :: id
        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)
    end subroutine end_id

end program test_zmumps_jobs_extra
