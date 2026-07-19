! Parallel analysis coverage for ZMUMPS — sequential Scotch vs PT-Scotch.
! Complex Hermitian-positive-definite mirror of
! test_dmumps_ptscotch_analysis: the same distributed sparse grid graph
! (ICNTL(18)=3), solved once through sequential Scotch (ICNTL(28)=1,
! ICNTL(7)=3) and once through PT-Scotch parallel analysis (ICNTL(28)=2,
! ICNTL(29)=1). See the d-variant header for the full rationale.

program test_zmumps_ptscotch_analysis
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_sparse_hpd_problem_z
    use target_mumps,          only: target_name, zmumps_struc, q2t_c
    use mumps_lifecycle,       only: mumps_begin, mumps_solve6, &
                                     mumps_check_solution, mumps_end
    use mpi
    implicit none

    integer, parameter :: nx = 8, ny = 8
    integer            :: n, ierr, nz, myid, nprocs, lo, hi, nz_loc
    logical            :: is_host
    complex(ep), allocatable :: x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    is_host = (myid == 0)
    call report_init('test_zmumps_ptscotch_analysis', target_name)

    call gen_sparse_hpd_problem_z(nx, ny, n, x_true, b, irn, jcn, A_trip, nz, &
                                  seed = 27101)

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

    subroutine run_case(label, parallel_ana)
        character(len=*), intent(in) :: label
        logical,          intent(in) :: parallel_ana

        call mumps_begin(id, MPI_COMM_WORLD, 0, ' ('//label//')')
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
        allocate(id%A_loc(nz_loc));    id%A_loc   = q2t_c(A_trip(lo:hi))
        if (is_host) then
            allocate(id%RHS(n));   id%RHS = q2t_c(b)
        end if

        call mumps_solve6(id, ' ('//label//')')

        if (is_host) call mumps_check_solution(id, x_true, label)

        call mumps_end(id)
    end subroutine run_case

end program test_zmumps_ptscotch_analysis
