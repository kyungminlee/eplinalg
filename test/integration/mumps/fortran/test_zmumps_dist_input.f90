! ICNTL(18)=3 — user-supplied distributed assembled input for ZMUMPS.
! Complex mirror of test_dmumps_dist_input: every rank contributes a
! disjoint contiguous slice of the shared deterministic triplet list, so
! MUMPS assembles the full matrix from the distributed pieces. See the
! d-variant header for the full rationale.

program test_zmumps_dist_input
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, &
                                     zmumps_struc, target_xmumps, q2t_c
    use mumps_lifecycle,       only: mumps_begin, mumps_check_solution, &
                                     mumps_end
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz, myid, nprocs, lo, hi, nz_loc
    logical            :: is_host
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    is_host = (myid == 0)
    call report_init('test_zmumps_dist_input', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 25001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)

    lo     = (myid * nz) / nprocs + 1
    hi     = ((myid + 1) * nz) / nprocs
    nz_loc = hi - lo + 1
    if (nz_loc < 0) nz_loc = 0

    call mumps_begin(id, MPI_COMM_WORLD, 0)
    id%ICNTL(18) = 3

    id%N        = n
    id%NZ_loc   = nz_loc
    id%NNZ_loc  = int(nz_loc, kind=8)
    allocate(id%IRN_loc(nz_loc));  id%IRN_loc = irn(lo:hi)
    allocate(id%JCN_loc(nz_loc));  id%JCN_loc = jcn(lo:hi)
    allocate(id%A_loc(nz_loc));    id%A_loc   = q2t_c(A_trip(lo:hi))

    if (is_host) then
        allocate(id%RHS(n));       id%RHS = q2t_c(b)
    end if

    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(18)=3 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if

    if (is_host) call mumps_check_solution(id, x_true, 'icntl18=3')

    call mumps_end(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_dist_input
