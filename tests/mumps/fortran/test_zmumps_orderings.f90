! ICNTL(7) ordering coverage for ZMUMPS — mirror of test_dmumps_orderings.

program test_zmumps_orderings
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_sparse_hpd_problem_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: nx = 8, ny = 8
    integer, parameter :: orderings(*) = [0, 2, 3, 4, 5, 6, 7]
    integer            :: n, ierr, i, nz, ord
    complex(ep), allocatable :: x_true(:), b(:), x_solve(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_orderings', target_name)

    call gen_sparse_hpd_problem_z(nx, ny, n, x_true, b, irn, jcn, A_trip, nz, &
                                  seed = 24001)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    do i = 1, size(orderings)
        ord = orderings(i)

        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_xmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
        id%ICNTL(7) = ord

        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));   id%RHS = q2t_c(b)

        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'JOB=6 (icntl7=', ord, &
                ') failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        write(label, '(a,i0)') 'icntl7=', ord
        call report_case(trim(label), err, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS, x_solve)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)
    end do

    deallocate(x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_orderings
