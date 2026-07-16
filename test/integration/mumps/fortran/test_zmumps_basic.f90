! Basic smoke + tolerance baseline for the migrated ZMUMPS via the
! kind16 ${LIB_PREFIX}mumps archive (called as `xmumps` here through
! the target_mumps wrapper). Mirrors test_dmumps_basic.f90 but with
! complex-valued coefficient and right-hand side.

program test_zmumps_basic
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: ns(*) = [8, 32]
    integer            :: ierr, i, n, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_basic', target_name)

    do i = 1, size(ns)
        n = ns(i)
        call gen_dense_problem_z(n, A, x_true, b, seed = 7001 + 31 * i)
        call dense_to_triplet_z (A, irn, jcn, A_trip, nz)

        id%COMM = MPI_COMM_WORLD
        id%PAR  = 1
        id%SYM  = 0
        id%JOB  = -1
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0)') 'JOB=-1 failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0

        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));    id%IRN = irn
        allocate(id%JCN(nz));    id%JCN = jcn
        allocate(id%A(nz));      id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));     id%RHS = q2t_c(b)

        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'JOB=6 failed, INFOG(1)=', &
                id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if

        allocate(x_solve(n))
        x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)

        deallocate(A, x_true, b, x_solve, irn, jcn, A_trip)
    end do

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_basic
