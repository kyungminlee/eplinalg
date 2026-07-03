! ICNTL(18)=3 — user-supplied distributed assembled matrix entry.
!
! With ICNTL(18)=3 each rank hands MUMPS its OWN share of the matrix via
! NZ_loc / NNZ_loc / IRN_loc / JCN_loc / A_loc, and MUMPS assembles the
! global matrix from the union of the per-rank contributions (duplicate
! (i,j) entries are summed, exactly like the centralized triplet path).
!
! This driver genuinely exercises that distributed assembly under
! ``mpiexec -n N``: every rank regenerates the full triplet list from the
! shared deterministic seed, then contributes a DISJOINT contiguous slice
! of it. The union over ranks is the full matrix with no duplication, so
! the assembled system is identical to the centralized one and the
! recovered solution must match x_true. At np=1 the single rank owns the
! whole slice, so this collapses to the original single-rank behaviour.
!
! N is set on every rank; the analysis/ordering still runs centrally
! (ICNTL(18)=3 changes only how the entries are supplied at JOB=1). The
! RHS and the centralized solution go through id%RHS on the host only.

program test_dmumps_dist_input
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps, &
                                     q2t_r, t2q_r
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz, myid, nprocs, lo, hi, nz_loc
    logical            :: is_host
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_solve(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    is_host = (myid == 0)
    call report_init('test_dmumps_dist_input', target_name)

    ! Every rank builds the full triplet (deterministic seed) and then
    ! owns the contiguous slice [lo:hi]. The balanced split gives every
    ! entry to exactly one rank, so the distributed union has no duplicates.
    call gen_dense_problem(n, A, x_true, b, seed = 14001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    lo     = (myid * nz) / nprocs + 1
    hi     = ((myid + 1) * nz) / nprocs
    nz_loc = hi - lo + 1
    if (nz_loc < 0) nz_loc = 0

    id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
    call target_qmumps(id)
    id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    id%ICNTL(18) = 3

    id%N        = n
    id%NZ_loc   = nz_loc
    id%NNZ_loc  = int(nz_loc, kind=8)
    allocate(id%IRN_loc(nz_loc));  id%IRN_loc = irn(lo:hi)
    allocate(id%JCN_loc(nz_loc));  id%JCN_loc = jcn(lo:hi)
    allocate(id%A_loc(nz_loc));    id%A_loc   = q2t_r(A_trip(lo:hi))

    ! RHS still goes through the centralized id%RHS path on the host.
    if (is_host) then
        allocate(id%RHS(n));   id%RHS = q2t_r(b)
    end if

    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0,a,i0)') 'JOB=6 with ICNTL(18)=3 failed, INFOG(1)=', &
            id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
        error stop 1
    end if

    if (is_host) then
        allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        call report_case('icntl18=3', err, tol)
        deallocate(id%RHS, x_solve)
        nullify(id%RHS)
    end if

    deallocate(id%IRN_loc, id%JCN_loc, id%A_loc)
    nullify(id%IRN_loc, id%JCN_loc, id%A_loc)
    id%JOB = -2;  call target_qmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_dist_input
