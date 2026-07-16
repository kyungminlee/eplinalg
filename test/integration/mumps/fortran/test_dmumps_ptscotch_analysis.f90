! Parallel analysis coverage for DMUMPS — sequential Scotch vs PT-Scotch.
!
! MUMPS offers two analysis (ordering) modes selected by ICNTL(28):
!   ICNTL(28)=1 — SEQUENTIAL analysis on the host; the ordering library is
!                 chosen by ICNTL(7) (here 3 = Scotch, the vendored
!                 sequential libscotch).
!   ICNTL(28)=2 — PARALLEL analysis distributed across the communicator;
!                 the ordering library is chosen by ICNTL(29) (here 1 =
!                 PT-Scotch, the distributed libptscotch increment). This
!                 is the ICNTL(28)=2 + ICNTL(29)=1 path that libptscotch
!                 exists to serve — inactive until the distributed archive
!                 is linked (before that MUMPS returns INFO(1)=-38).
!
! Both passes solve the SAME sparse SPD system supplied as a DISTRIBUTED
! assembled matrix (ICNTL(18)=3): every rank contributes a disjoint slice
! of the 2D-Laplacian triplet list, MUMPS assembles the union. The matrix
! is the sparse grid graph gen_sparse_spd_problem builds — nested
! dissection (sequential Scotch AND PT-Scotch alike) needs real vertex
! separators, which a dense graph would not provide. Each pass must
! recover x_true to the quad-precision tolerance, proving both orderings
! drive the factorization to the same solution.
!
! Runs under ``mpiexec -n N``; PT-Scotch genuinely distributes the graph
! across the ranks. At np=1 both passes collapse to a single-rank solve.

program test_dmumps_ptscotch_analysis
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_sparse_spd_problem
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps, &
                                     q2t_r, t2q_r
    use mpi
    implicit none

    integer, parameter :: nx = 8, ny = 8
    integer            :: n, ierr, nz, myid, nprocs, lo, hi, nz_loc
    logical            :: is_host
    real(ep), allocatable :: x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc)    :: id
    real(ep)              :: tol

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    is_host = (myid == 0)
    call report_init('test_dmumps_ptscotch_analysis', target_name)

    ! Every rank builds the full sparse triplet from the shared seed, then
    ! owns the contiguous slice [lo:hi]; the union over ranks is the full
    ! matrix with no duplicated entries.
    call gen_sparse_spd_problem(nx, ny, n, x_true, b, irn, jcn, A_trip, nz, &
                                seed = 7101)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    lo     = (myid * nz) / nprocs + 1
    hi     = ((myid + 1) * nz) / nprocs
    nz_loc = hi - lo + 1
    if (nz_loc < 0) nz_loc = 0

    ! Pass A: sequential analysis, Scotch (ICNTL(28)=1, ICNTL(7)=3).
    call run_case('dist+scotch',   parallel_ana = .false.)
    ! Pass B: parallel analysis, PT-Scotch (ICNTL(28)=2, ICNTL(29)=1).
    call run_case('dist+ptscotch', parallel_ana = .true.)

    deallocate(x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    ! One full JOB=-1 .. 6 .. -2 solve of the distributed system with the
    ! selected analysis mode, reporting the host-side relative error.
    subroutine run_case(label, parallel_ana)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: parallel_ana
        real(ep), allocatable :: x_solve(:)
        real(ep)              :: err

        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0)') 'JOB=-1 (', label, ') failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
        id%ICNTL(18) = 3                      ! distributed assembled entry
        if (parallel_ana) then
            id%ICNTL(28) = 2                  ! parallel analysis
            id%ICNTL(29) = 1                  ! ... via PT-Scotch
        else
            id%ICNTL(28) = 1                  ! sequential analysis
            id%ICNTL(7)  = 3                  ! ... via Scotch
        end if

        id%N       = n
        id%NZ_loc  = nz_loc
        id%NNZ_loc = int(nz_loc, kind=8)
        allocate(id%IRN_loc(nz_loc));  id%IRN_loc = irn(lo:hi)
        allocate(id%JCN_loc(nz_loc));  id%JCN_loc = jcn(lo:hi)
        allocate(id%A_loc(nz_loc));    id%A_loc   = q2t_r(A_trip(lo:hi))
        if (is_host) then
            allocate(id%RHS(n));   id%RHS = q2t_r(b)
        end if

        id%JOB = 6
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0,a,i0)') 'JOB=6 (', label, ') failed, INFOG(1)=', &
                id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if

        if (is_host) then
            allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
            err = max_rel_err_vec(x_solve, x_true)
            call report_case(label, err, tol)
            deallocate(id%RHS, x_solve)
            nullify(id%RHS)
        end if

        deallocate(id%IRN_loc, id%JCN_loc, id%A_loc)
        nullify(id%IRN_loc, id%JCN_loc, id%A_loc)
        id%JOB = -2;  call target_qmumps(id)
    end subroutine run_case

end program test_dmumps_ptscotch_analysis
