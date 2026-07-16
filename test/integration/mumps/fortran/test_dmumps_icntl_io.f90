! ICNTL coverage for input/output formats and scaling — exercises
! struct fields the migrated archive renamed but no other test
! touches:
!
!   ICNTL(8) — scaling strategy. Values cover the SCALING dispatcher
!     (qmumps_driver.F:1284-1623) which fills COLSCA / ROWSCA arrays
!     (REAL(KIND=16) on the migrated side).
!   ICNTL(20) — RHS format. ICNTL(20)=1 means sparse RHS, requiring
!     RHS_SPARSE / IRHS_SPARSE / IRHS_PTR / NZ_RHS instead of dense
!     RHS. Different code path through the solve driver.
!   ICNTL(21) — solution distribution. ICNTL(21)=1 returns the
!     solution distributed via SOL_loc / ISOL_loc / LSOL_loc on each
!     rank, INFO(23) components apiece. At np=1 the lone rank owns all
!     N; at np≥2 the slice is genuinely partitioned across ranks.

program test_dmumps_icntl_io
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
    integer            :: ierr, nz, myid
    integer, parameter :: scaling_modes(*) = [0, 1, 7, 77]  ! none / diag / sym / auto
    integer            :: i, sc
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    real(ep), allocatable :: x_solve(:)
    type(dmumps_struc)    :: id
    real(ep)              :: err, tol
    character(len=48)     :: label

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call report_init('test_dmumps_icntl_io', target_name)
    call gen_dense_problem(n, A, x_true, b, seed = 12001)
    call dense_to_triplet (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    ! ── ICNTL(8) scaling strategies ────────────────────────────────
    do i = 1, size(scaling_modes)
        sc = scaling_modes(i)
        call init_id(id)
        id%ICNTL(8) = sc
        call attach_dense(id, n, nz, irn, jcn, A_trip, b)
        id%JOB = 6
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'ICNTL(8)=', sc, &
                ' failed: INFOG(1)=', id%INFOG(1)
            error stop 1
        end if
        allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        write(label, '(a,i0)') 'icntl8=', sc
        call report_case(trim(label), err, tol)
        call end_id(id);  deallocate(x_solve)
    end do

    ! ── ICNTL(21)=1 distributed solution ───────────────────────────
    !    Each rank must allocate SOL_loc / ISOL_loc to LSOL_loc (set to
    !    N here, an upper bound on any rank's local slice). At np=1 the
    !    lone rank owns all N components; at np≥2 the slice is genuinely
    !    partitioned across ranks (see the INFO(23) note below).
    call init_id(id)
    id%ICNTL(21) = 1
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    ! For ICNTL(21)=1 we need SOL_loc / ISOL_loc allocated; pre-set
    ! LSOL_loc so MUMPS knows the user-side capacity.
    id%LSOL_loc = n
    allocate(id%SOL_loc(n))
    allocate(id%ISOL_loc(n))
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(21)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    ! Distributed solution (ICNTL(21)=1): each rank receives its own
    ! INFO(23) components in SOL_loc, with the global row indices in
    ! ISOL_loc; the union over ranks is the full solution. INFO(23) (the
    ! number of pivots eliminated on this processor, per the User Guide)
    ! is the documented size of the local slice — the struct's NSOL_loc
    ! field is NOT written by the solve. At np=1 the single rank owns all
    ! n components, so this collapses to the centralized check.
    ! max_rel_err_vec is max|Δ|/max|x_true|: the denominator is global
    ! (x_true is full on every rank) and the numerator is the max over
    ! ranks of each rank's local max|Δ|, so a single scalar MPI_MAX
    ! suffices. It is reduced in double precision — the error magnitude
    ! (~1e-19) is well inside double's exponent range and the tolerance
    ! has ~15 significant digits of head-room, so no extended-precision
    ! MPI datatype is needed.
    block
        real(ep) :: local_num, denom
        real(8)  :: loc8, glob8
        integer  :: idx, gi
        local_num = 0.0_ep
        do idx = 1, id%INFO(23)
            gi = id%ISOL_loc(idx)
            local_num = max(local_num, &
                            abs(t2q_r(id%SOL_loc(idx)) - x_true(gi)))
        end do
        denom = maxval(abs(x_true))
        loc8  = real(local_num, kind=8)
        glob8 = 0.0d0
        call MPI_Reduce(loc8, glob8, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                        0, MPI_COMM_WORLD, ierr)
        if (denom < tiny(1.0_ep)) then
            err = real(glob8, kind=ep)
        else
            err = real(glob8, kind=ep) / denom
        end if
    end block
    call report_case('icntl21=1', err, tol)
    deallocate(id%SOL_loc, id%ISOL_loc)
    call end_id(id)

    ! ── ICNTL(20)=1 sparse RHS — single-column sparse RHS that
    !    represents the same dense RHS we'd otherwise pass. NZ_RHS
    !    is the total nonzero count across all RHS columns. For
    !    NRHS=1 the per-column structure is degenerate (just
    !    IRHS_PTR=[1, NZ_RHS+1]), but exercising the dispatch is
    !    what we want.
    !
    !    Centralized sparse RHS is a HOST-ONLY input: the solve driver
    !    (?sol_driver.F) rejects with INFO(1)=-22 any non-master rank
    !    whose RHS / RHS_SPARSE / IRHS_SPARSE / IRHS_PTR are associated.
    !    So only the host allocates the sparse-RHS arrays and the dense
    !    solution vector RHS; the slaves leave them null. The matrix
    !    (IRN/JCN/A) stays centralized on every rank as elsewhere. ───
    call init_id(id)
    id%ICNTL(20) = 1
    id%N    = n
    id%NNZ  = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = q2t_r(A_trip)
    id%NRHS    = 1
    if (myid == 0) then
        ! Sparse RHS layout for NRHS=1: every entry of b is "non-zero"
        ! (we don't actually sparsify it — that defeats the test) and
        ! IRHS_SPARSE just enumerates 1..n, IRHS_PTR is [1, n+1].
        id%LRHS    = n
        id%NZ_RHS  = n
        allocate(id%RHS_SPARSE(n));   id%RHS_SPARSE  = q2t_r(b)
        allocate(id%IRHS_SPARSE(n))
        allocate(id%IRHS_PTR(2))
        do i = 1, n;  id%IRHS_SPARSE(i) = i;  end do
        id%IRHS_PTR(1) = 1
        id%IRHS_PTR(2) = n + 1
        ! Solution still lands in RHS (centralized) — host only.
        allocate(id%RHS(n));  id%RHS = q2t_r(0.0_ep)
    end if
    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(20)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    if (myid == 0) then
        allocate(x_solve(n));  x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        deallocate(id%RHS_SPARSE, id%IRHS_SPARSE, id%IRHS_PTR, x_solve)
    end if
    call report_case('icntl20=1', err, tol)
    call end_id(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine init_id(id)
        type(dmumps_struc), intent(inout) :: id
        id%COMM = MPI_COMM_WORLD
        id%PAR  = 1
        id%SYM  = 0
        id%JOB  = -1
        call target_qmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    end subroutine init_id

    subroutine attach_dense(id, n, nz, irn, jcn, A_trip, b)
        type(dmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        real(ep),           intent(in)    :: A_trip(:), b(:)
        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = q2t_r(A_trip)
        allocate(id%RHS(n));   id%RHS = q2t_r(b)
    end subroutine attach_dense

    subroutine end_id(id)
        type(dmumps_struc), intent(inout) :: id
        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        id%JOB = -2
        call target_qmumps(id)
    end subroutine end_id

end program test_dmumps_icntl_io
