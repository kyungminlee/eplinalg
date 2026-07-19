! Shared MUMPS instance-lifecycle helpers for the hand-written test
! drivers under fortran/. Every solve-check driver runs the same
! skeleton — JOB=-1 init + INFOG abort, ICNTL(1..4) silencing,
! centralized triplet/RHS injection, JOB=6 + INFOG abort, host-side
! solution check, deallocate/nullify, JOB=-2 — so it lives here once,
! per struc type via generic dispatch on dmumps_struc / zmumps_struc.
! Drivers keep only their distinctive setup (problem generator, extra
! ICNTL overrides, case loop, custom checks).
!
! MPI is deliberately NOT used here: drivers pass their communicator
! (MPI_COMM_WORLD) into mumps_begin, so this module compiles without
! mpi.mod and stays neutral across the plain/_seq link variants.
!
! `ctx` is an optional pre-formatted context fragment (e.g.
! ' (sym=1)') spliced into the abort diagnostics; the abort itself is
! collective-consistent because MUMPS broadcasts INFOG to all ranks.
module mumps_lifecycle
    use prec_kinds,   only: ep
    use compare,      only: max_rel_err_vec, max_rel_err_vec_z
    use prec_report,  only: report_case
    use target_mumps, only: target_eps, &
                            dmumps_struc, target_qmumps, q2t_r, t2q_r, &
                            zmumps_struc, target_xmumps, q2t_c, t2q_c
    implicit none
    private
    public :: mumps_begin, mumps_load_triplet, mumps_solve6, &
              mumps_check_solution, mumps_end, mumps_default_tol

    ! JOB=-1 init (COMM/PAR/SYM) + INFOG abort + diagnostics silencing.
    interface mumps_begin
        module procedure mumps_begin_r, mumps_begin_c
    end interface mumps_begin

    ! Centralized assembled input: N/NNZ + IRN/JCN/A/RHS injection.
    interface mumps_load_triplet
        module procedure mumps_load_triplet_r, mumps_load_triplet_c
    end interface mumps_load_triplet

    ! Combined analysis + factorization + solve (JOB=6) + INFOG abort.
    interface mumps_solve6
        module procedure mumps_solve6_r, mumps_solve6_c
    end interface mumps_solve6

    ! Host-side check of the centralized solution in id%RHS against
    ! x_true; only call where id%RHS is populated (the host rank).
    interface mumps_check_solution
        module procedure mumps_check_solution_r, mumps_check_solution_c
    end interface mumps_check_solution

    ! Guarded deallocate/nullify of the user-owned pointer arrays
    ! (centralized and _loc) + JOB=-2. Safe on ranks that never
    ! allocated (JOB=-1 nullifies every pointer component).
    interface mumps_end
        module procedure mumps_end_r, mumps_end_c
    end interface mumps_end

contains

    ! Repo-standard n**3-scaled solve tolerance (matches
    ! tests/lapack/linear_solve/test_dgesv.f90).
    pure function mumps_default_tol(n) result(tol)
        integer, intent(in) :: n
        real(ep) :: tol
        tol = 16.0_ep * real(n, ep)**3 * target_eps
    end function mumps_default_tol

    subroutine mumps_begin_r(id, comm, sym, ctx)
        type(dmumps_struc), intent(inout)        :: id
        integer,            intent(in)           :: comm, sym
        character(len=*),   intent(in), optional :: ctx

        id%COMM = comm
        id%PAR  = 1
        id%SYM  = sym
        id%JOB  = -1
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0)') 'JOB=-1', ctx_str(ctx), &
                ' failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        ! Silence MUMPS diagnostic / info output.
        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0
    end subroutine mumps_begin_r

    subroutine mumps_begin_c(id, comm, sym, ctx)
        type(zmumps_struc), intent(inout)        :: id
        integer,            intent(in)           :: comm, sym
        character(len=*),   intent(in), optional :: ctx

        id%COMM = comm
        id%PAR  = 1
        id%SYM  = sym
        id%JOB  = -1
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0)') 'JOB=-1', ctx_str(ctx), &
                ' failed, INFOG(1)=', id%INFOG(1)
            error stop 1
        end if

        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0
    end subroutine mumps_begin_c

    subroutine mumps_load_triplet_r(id, n, nz, irn, jcn, a_trip, b)
        type(dmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        real(ep),           intent(in)    :: a_trip(:), b(:)

        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));    id%IRN = irn
        allocate(id%JCN(nz));    id%JCN = jcn
        allocate(id%A(nz));      id%A   = q2t_r(a_trip)
        allocate(id%RHS(n));     id%RHS = q2t_r(b)
    end subroutine mumps_load_triplet_r

    subroutine mumps_load_triplet_c(id, n, nz, irn, jcn, a_trip, b)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        complex(ep),        intent(in)    :: a_trip(:), b(:)

        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));    id%IRN = irn
        allocate(id%JCN(nz));    id%JCN = jcn
        allocate(id%A(nz));      id%A   = q2t_c(a_trip)
        allocate(id%RHS(n));     id%RHS = q2t_c(b)
    end subroutine mumps_load_triplet_c

    subroutine mumps_solve6_r(id, ctx)
        type(dmumps_struc), intent(inout)        :: id
        character(len=*),   intent(in), optional :: ctx

        id%JOB = 6
        call target_qmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0,a,i0)') 'JOB=6', ctx_str(ctx), &
                ' failed, INFOG(1)=', id%INFOG(1), &
                ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if
    end subroutine mumps_solve6_r

    subroutine mumps_solve6_c(id, ctx)
        type(zmumps_struc), intent(inout)        :: id
        character(len=*),   intent(in), optional :: ctx

        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(3a,i0,a,i0)') 'JOB=6', ctx_str(ctx), &
                ' failed, INFOG(1)=', id%INFOG(1), &
                ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if
    end subroutine mumps_solve6_c

    subroutine mumps_check_solution_r(id, x_true, label, tol)
        type(dmumps_struc), intent(in)           :: id
        real(ep),           intent(in)           :: x_true(:)
        character(len=*),   intent(in)           :: label
        real(ep),           intent(in), optional :: tol
        real(ep) :: x_solve(size(x_true)), err, t

        x_solve = t2q_r(id%RHS)
        err = max_rel_err_vec(x_solve, x_true)
        t = mumps_default_tol(size(x_true))
        if (present(tol)) t = tol
        call report_case(label, err, t)
    end subroutine mumps_check_solution_r

    subroutine mumps_check_solution_c(id, x_true, label, tol)
        type(zmumps_struc), intent(in)           :: id
        complex(ep),        intent(in)           :: x_true(:)
        character(len=*),   intent(in)           :: label
        real(ep),           intent(in), optional :: tol
        complex(ep) :: x_solve(size(x_true))
        real(ep)    :: err, t

        x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        t = mumps_default_tol(size(x_true))
        if (present(tol)) t = tol
        call report_case(label, err, t)
    end subroutine mumps_check_solution_c

    subroutine mumps_end_r(id)
        type(dmumps_struc), intent(inout) :: id

        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        ! Defensive: JOB=-1 nullifies these too, but keep the
        ! disassociation explicit rather than relying on it.
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        if (associated(id%IRN_loc)) deallocate(id%IRN_loc)
        if (associated(id%JCN_loc)) deallocate(id%JCN_loc)
        if (associated(id%A_loc))   deallocate(id%A_loc)
        nullify(id%IRN_loc, id%JCN_loc, id%A_loc)

        id%JOB = -2
        call target_qmumps(id)
    end subroutine mumps_end_r

    subroutine mumps_end_c(id)
        type(zmumps_struc), intent(inout) :: id

        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        if (associated(id%IRN_loc)) deallocate(id%IRN_loc)
        if (associated(id%JCN_loc)) deallocate(id%JCN_loc)
        if (associated(id%A_loc))   deallocate(id%A_loc)
        nullify(id%IRN_loc, id%JCN_loc, id%A_loc)

        id%JOB = -2
        call target_xmumps(id)
    end subroutine mumps_end_c

    function ctx_str(ctx) result(s)
        character(len=*), intent(in), optional :: ctx
        character(len=:), allocatable :: s
        if (present(ctx)) then
            s = ctx
        else
            s = ''
        end if
    end function ctx_str

end module mumps_lifecycle
