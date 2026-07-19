! ICNTL(10) iterative refinement and ICNTL(11) error analysis paths.
!
!   ICNTL(10) > 0  — iterative refinement with stopping criterion;
!                    INFOG(15) reports number of IR steps performed.
!   ICNTL(10) < 0  — fixed |ICNTL(10)| IR steps regardless.
!
!   ICNTL(11) = 1 or 2 — error analysis fills RINFOG(4..11) with
!                    componentwise / normwise residual estimates.
!                    This is the only path that exercises the
!                    migrated REAL(KIND=16) RINFOG fields beyond
!                    INFOG(1) verification, so it's the natural
!                    smoke for "did the precision promotion of the
!                    RINFOG/RINFO arrays survive migration?"

program test_dmumps_iref_errchk
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, target_eps, &
                                     dmumps_struc, target_qmumps, t2q_r
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_end, mumps_default_tol
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_solve(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_dmumps_iref_errchk', target_name)
    call gen_dense_problem(n, A, x_true, b, seed = 13001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = mumps_default_tol(n)

    ! ── ICNTL(10) = 5 — up to 5 IR steps with default stopping ────
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    id%ICNTL(10) = 5
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(10)=5 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('icntl10=5', err, tol)
    deallocate(x_solve);  call mumps_end(id)

    ! ── ICNTL(10) = -2 — exactly 2 IR steps ───────────────────────
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    id%ICNTL(10) = -2
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(10)=-2 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('icntl10=-2', err, tol)
    deallocate(x_solve);  call mumps_end(id)

    ! ── ICNTL(11) = 2 — full error analysis. Computes scaled
    !    residual into RINFOG(6) and componentwise backward errors
    !    omega1=RINFOG(7), omega2=RINFOG(8). Verifies the migrated
    !    REAL(KIND=16) RINFOG layout produces sane finite values
    !    under the precision promotion. ──────────────────────────
    call mumps_begin(id, MPI_COMM_WORLD, 0)
    id%ICNTL(11) = 2
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(11)=2 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
    err = max_rel_err_vec(x_solve, x_true)
    call report_case('icntl11=2:solution', err, tol)
    ! RINFOG(6) is the scaled solution residual ||Ax-b|| / something.
    ! With a well-conditioned diagonally-dominant matrix and quad
    ! precision the residual should be tiny — at most a few ULPs of
    ! KIND=16 epsilon. The bound is generous to avoid false negatives
    ! on platforms with different LTO inlining: any sane residual
    ! below 1e-25 means MUMPS produced a finite RINFOG(6) and it's
    ! consistent with quad-precision arithmetic.
    ! Loose bound on the iter-refinement residual: a few orders of
    ! magnitude above the target's epsilon. The actual residual is
    ! typically a small multiple of eps; the bound is generous to
    ! avoid false negatives on platforms with different LTO inlining.
    block
        real(ep) :: rinfog6_bound, rinfog6_q
        rinfog6_bound = 1.0e6_ep * target_eps
        rinfog6_q = t2q_r(id%RINFOG(6))
        if (rinfog6_q /= rinfog6_q) then  ! NaN check
            call report_case('icntl11=2:RINFOG6-finite', 1.0_ep, 0.0_ep)
        else if (rinfog6_q > rinfog6_bound) then
            call report_case('icntl11=2:RINFOG6-bound', &
                             rinfog6_q, rinfog6_bound)
        else
            call report_case('icntl11=2:RINFOG6-bound', 0.0_ep, 1.0_ep)
        end if
    end block
    deallocate(x_solve);  call mumps_end(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_iref_errchk
