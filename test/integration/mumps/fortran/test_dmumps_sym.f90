! Symmetry-flag coverage for DMUMPS:
!   SYM = 0 — general unsymmetric (LU)
!   SYM = 1 — symmetric positive-definite (Cholesky)
!   SYM = 2 — general symmetric (LDL^T, Bunch-Kaufman)
!
! For each SYM value, generate a quad-precision dense problem with
! the matching structure, hand the upper triangle (SYM > 0) or the
! full matrix (SYM = 0) to MUMPS in triplet form, and check the
! recovered solution against x_true at REAL(KIND=ep).

program test_dmumps_sym
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, gen_spd_dense_problem, &
                                     gen_general_sym_problem, &
                                     dense_to_triplet, dense_to_sym_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps, &
                                     q2t_r, t2q_r
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, sym, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol
    character(len=48)     :: label

    call MPI_INIT(ierr)
    call report_init('test_dmumps_sym', target_name)

    do sym = 0, 2
        select case (sym)
        case (0); call gen_dense_problem        (n, A, x_true, b, seed = 2001)
        case (1); call gen_spd_dense_problem    (n, A, x_true, b, seed = 2002)
        case (2); call gen_general_sym_problem  (n, A, x_true, b, seed = 2003)
        end select

        if (sym == 0) then
            call dense_to_triplet    (A, irn, jcn, A_trip, nz)
        else
            call dense_to_sym_triplet(A, irn, jcn, A_trip, nz)
        end if

        id%COMM = MPI_COMM_WORLD
        id%PAR  = 1
        id%SYM  = sym
        id%JOB  = -1
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'JOB=-1 (sym=', sym, ') failed, INFOG(1)=', id%INFOG(1)
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
        allocate(id%A(nz));      id%A   = q2t_r(A_trip)
        allocate(id%RHS(n));     id%RHS = q2t_r(b)

        id%JOB = 6
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0,a,i0)') 'JOB=6 (sym=', sym, ') failed, INFOG(1)=', &
                id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if

        allocate(x_solve(n))
        x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'sym=', sym
        call report_case(trim(label), err, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_qmumps(id)

        deallocate(A, x_true, b, x_solve, irn, jcn, A_trip)
    end do

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_sym
