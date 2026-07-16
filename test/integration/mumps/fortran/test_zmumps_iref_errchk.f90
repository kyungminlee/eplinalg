! ICNTL(10)/(11) coverage for ZMUMPS — mirror of test_dmumps_iref_errchk.

program test_zmumps_iref_errchk
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c, t2q_r
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_solve(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_zmumps_iref_errchk', target_name)
    call gen_dense_problem_z(n, A, x_true, b, seed = 27001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    call init_id(id)
    id%ICNTL(10) = 5
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_xmumps(id)
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('icntl10=5', err, tol)
    deallocate(x_solve);  call end_id(id)

    call init_id(id)
    id%ICNTL(10) = -2
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_xmumps(id)
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('icntl10=-2', err, tol)
    deallocate(x_solve);  call end_id(id)

    call init_id(id)
    id%ICNTL(11) = 2
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_xmumps(id)
    allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
    err = max_rel_err_vec_z(x_solve, x_true)
    call report_case('icntl11=2:solution', err, tol)
    ! See test_dmumps_iref_errchk for the rationale on this bound.
    block
        real(ep) :: rinfog6_bound, rinfog6_q
        rinfog6_bound = 1.0e6_ep * target_eps
        rinfog6_q = t2q_r(id%RINFOG(6))
        if (rinfog6_q /= rinfog6_q) then
            call report_case('icntl11=2:RINFOG6-finite', 1.0_ep, 0.0_ep)
        else if (rinfog6_q > rinfog6_bound) then
            call report_case('icntl11=2:RINFOG6-bound', &
                             rinfog6_q, rinfog6_bound)
        else
            call report_case('icntl11=2:RINFOG6-bound', 0.0_ep, 1.0_ep)
        end if
    end block
    deallocate(x_solve);  call end_id(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine init_id(id)
        type(zmumps_struc), intent(inout) :: id
        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_xmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    end subroutine init_id

    subroutine attach_dense(id, n, nz, irn, jcn, A_trip, b)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        complex(ep),        intent(in)    :: A_trip(:), b(:)
        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));   id%RHS = q2t_c(b)
    end subroutine attach_dense

    subroutine end_id(id)
        type(zmumps_struc), intent(inout) :: id
        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2;  call target_xmumps(id)
    end subroutine end_id

end program test_zmumps_iref_errchk
