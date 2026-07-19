! JOB-phasing equivalence for DMUMPS:
!   Combined  JOB=6 (analyze + factor + solve)  in one call
!   Phased    JOB=-1, 1, 2, 3, -2 in five calls
! The two paths must produce the same numerical solution to within
! the per-case tolerance — INFOG(1) must be 0 from each, and the
! RHS at exit must match.

program test_dmumps_jobs
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, &
                                     dmumps_struc, target_qmumps, t2q_r
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_end, mumps_default_tol
    use mpi
    implicit none

    integer, parameter :: n = 12
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_combined(:), x_phased(:)
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_dmumps_jobs', target_name)

    call gen_dense_problem(n, A, x_true, b, seed = 3001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)

    ! ── Path A: combined JOB=6 ──────────────────────────────────────
    call mumps_solve_jobs(jobs = [6], n=n, nz=nz, irn=irn, jcn=jcn, &
                          A_trip=A_trip, b=b, x_solve=x_combined)

    ! ── Path B: phased JOB=1 → 2 → 3 ────────────────────────────────
    call mumps_solve_jobs(jobs = [1, 2, 3], n=n, nz=nz, irn=irn, jcn=jcn, &
                          A_trip=A_trip, b=b, x_solve=x_phased)

    ! Each path against the quad ground-truth.
    err = max_rel_err_vec(x_combined, x_true)
    tol = mumps_default_tol(n)
    call report_case('combined-vs-truth', err, tol)

    err = max_rel_err_vec(x_phased, x_true)
    call report_case('phased-vs-truth', err, tol)

    ! And the two paths against each other. JOB=6 dispatches the same
    ! analysis/factor/solve internals as JOB=1+2+3 on the same struct.
    ! On native real kinds the two paths order operations identically and
    ! agree bit-for-bit, but that is not a contract: extended-precision
    ! (multifloats) arithmetic accumulates in a different order and the two
    ! paths can legitimately diverge by a few ULP. Gate on the same
    ! tolerance as the truth comparisons rather than demanding bit-equality.
    err = max_rel_err_vec(x_phased, x_combined)
    call report_case('phased-vs-combined', err, tol)

    deallocate(A, x_true, b, irn, jcn, A_trip, x_combined, x_phased)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    ! Run a JOB sequence (init = -1 prepended, end = -2 appended) and
    ! capture the solution after the LAST job that touched the RHS.
    subroutine mumps_solve_jobs(jobs, n, nz, irn, jcn, A_trip, b, x_solve)
        integer,  intent(in)               :: jobs(:), n, nz, irn(:), jcn(:)
        real(ep), intent(in)               :: A_trip(:), b(:)
        real(ep), allocatable, intent(out) :: x_solve(:)
        type(dmumps_struc) :: idl
        integer :: k

        call mumps_begin(idl, MPI_COMM_WORLD, 0)
        call mumps_load_triplet(idl, n, nz, irn, jcn, A_trip, b)

        do k = 1, size(jobs)
            idl%JOB = jobs(k)
            call target_qmumps(idl)
            if (idl%INFOG(1) < 0) then
                write(*, '(a,i0,a,i0)') 'JOB=', jobs(k), &
                    ' failed, INFOG(1)=', idl%INFOG(1)
                error stop 1
            end if
        end do

        allocate(x_solve(n))
        x_solve = t2q_r(idl%RHS)

        call mumps_end(idl)
    end subroutine mumps_solve_jobs

end program test_dmumps_jobs
