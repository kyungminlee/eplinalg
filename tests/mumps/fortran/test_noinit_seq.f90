! Sequential (libmpiseq) no-initialization solve check.
!
! Proves a seq-build consumer can drive the migrated MUMPS solver for
! BOTH supported types of this target (real + complex) with NEITHER
! MPI_Init NOR the per-target MPI-handle registration
! (multifloats_mpi_init / quad_mpi_init). This is the scenario a plain
! sequential application runs in: link libmpiseq, call ?mumps, never
! touch MPI.
!
! Historically this path STOPped inside libseq's MUMPS_COPY on
! multifloats: the MPI_FLOAT64X2 / MPI_COMPLEX64X2 datatype handles were
! 0 until multifloats_mpi_init ran, so the first non-in-place reduction
! hit an unknown datatype and libseq's MUMPS_COPY aborted. The datatype
! handles now default to the libmpiseq derived-type sentinel, so the
! solve completes with no init (runtime/multifloats-mpi/multifloats_mpi.cpp).
! kind16 / kind10 never depended on init for a single-rank solve (custom
! reduce ops are only reached at np >= 2), so they exercise the same
! no-init contract.
!
! Built ONLY in the _seq (libmpiseq) variant: the real-MPI variant
! legitimately requires MPI_Init, so a no-init driver is meaningless
! there. See the seq-only guard in tests/mumps/CMakeLists.txt.
!
! For each type a small dense system A x = b is solved (JOB=-1 -> 6 ->
! -2) and BOTH the solution accuracy and the residual norm are checked:
!   * solution   max|x_solve - x_true| / max|x_true|
!   * residual   max|b - A x_solve|    / max|b|
! against the repo-standard O(n^3) eps tolerance.

program test_noinit_seq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize, &
                               report_check_status
    use compare,         only: max_rel_err_vec, max_rel_err_vec_z
    use test_data_mumps, only: gen_dense_problem, dense_to_triplet, &
                               gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,    only: target_name, target_eps, &
                               dmumps_struc, zmumps_struc, &
                               target_qmumps_noinit, target_xmumps_noinit, &
                               q2t_r, t2q_r, q2t_c, t2q_c
    implicit none

    ! Small but non-trivial: an elimination tree with real merges, still
    ! cheap enough that host double precision generates the reference
    ! problem faithfully for every target.
    integer, parameter :: n = 8
    real(ep)           :: tol

    call report_init('test_noinit_seq', target_name)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    call solve_real()
    call solve_complex()

    call report_finalize()
    call report_check_status()

contains

    subroutine abort_if_failed(infog1, stage)
        integer,          intent(in) :: infog1
        character(len=*), intent(in) :: stage
        if (infog1 < 0) then
            write(*, '(a,a,a,i0)') 'MUMPS ', stage, ' failed, INFOG(1)=', infog1
            error stop 1
        end if
    end subroutine abort_if_failed

    subroutine solve_real()
        type(dmumps_struc)    :: id
        real(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:), r(:)
        integer,  allocatable :: irn(:), jcn(:)
        real(ep), allocatable :: A_trip(:)
        integer               :: nz
        real(ep)              :: err, resid

        call gen_dense_problem(n, A, x_true, b, seed = 7001)
        call dense_to_triplet(A, irn, jcn, A_trip, nz)

        ! Deliberately NO MPI_Init and NO multifloats_mpi_init /
        ! quad_mpi_init. A sequential consumer sets COMM to any value;
        ! libmpiseq ignores it entirely.
        id%COMM = 0
        id%PAR  = 1
        id%SYM  = 0
        id%JOB  = -1
        call target_qmumps_noinit(id)
        call abort_if_failed(id%INFOG(1), 'real JOB=-1')

        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0

        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz)); id%IRN = irn
        allocate(id%JCN(nz)); id%JCN = jcn
        allocate(id%A(nz));   id%A   = q2t_r(A_trip)
        allocate(id%RHS(n));  id%RHS = q2t_r(b)

        id%JOB = 6
        call target_qmumps_noinit(id)
        call abort_if_failed(id%INFOG(1), 'real JOB=6')

        allocate(x_solve(n)); x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        call report_case('real solution', err, tol)

        allocate(r(n)); r = b - matmul(A, x_solve)
        resid = maxval(abs(r)) / max(maxval(abs(b)), tiny(1.0_ep))
        call report_case('real residual', resid, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_qmumps_noinit(id)

        deallocate(A, x_true, b, x_solve, r, irn, jcn, A_trip)
    end subroutine solve_real

    subroutine solve_complex()
        type(zmumps_struc)       :: id
        complex(ep), allocatable :: A(:,:), x_true(:), b(:), x_solve(:), r(:)
        integer,     allocatable :: irn(:), jcn(:)
        complex(ep), allocatable :: A_trip(:)
        integer                  :: nz
        real(ep)                 :: err, resid

        call gen_dense_problem_z(n, A, x_true, b, seed = 7013)
        call dense_to_triplet_z(A, irn, jcn, A_trip, nz)

        id%COMM = 0
        id%PAR  = 1
        id%SYM  = 0
        id%JOB  = -1
        call target_xmumps_noinit(id)
        call abort_if_failed(id%INFOG(1), 'complex JOB=-1')

        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0

        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz)); id%IRN = irn
        allocate(id%JCN(nz)); id%JCN = jcn
        allocate(id%A(nz));   id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));  id%RHS = q2t_c(b)

        id%JOB = 6
        call target_xmumps_noinit(id)
        call abort_if_failed(id%INFOG(1), 'complex JOB=6')

        allocate(x_solve(n)); x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        call report_case('complex solution', err, tol)

        allocate(r(n)); r = b - matmul(A, x_solve)
        resid = maxval(abs(r)) / max(maxval(abs(b)), tiny(1.0_ep))
        call report_case('complex residual', resid, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps_noinit(id)

        deallocate(A, x_true, b, x_solve, r, irn, jcn, A_trip)
    end subroutine solve_complex

end program test_noinit_seq
