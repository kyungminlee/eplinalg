! ICNTL coverage for ZMUMPS — mirror of test_dmumps_icntl_io.

program test_zmumps_icntl_io
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz, myid
    integer, parameter :: scaling_modes(*) = [0, 1, 7, 77]
    integer            :: i, sc
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_solve(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call report_init('test_zmumps_icntl_io', target_name)
    call gen_dense_problem_z(n, A, x_true, b, seed = 26001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    do i = 1, size(scaling_modes)
        sc = scaling_modes(i)
        call init_id(id)
        id%ICNTL(8) = sc
        call attach_dense(id, n, nz, irn, jcn, A_trip, b)
        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'ICNTL(8)=', sc, &
                ' failed: INFOG(1)=', id%INFOG(1)
            error stop 1
        end if
        allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        write(label, '(a,i0)') 'icntl8=', sc
        call report_case(trim(label), err, tol)
        call end_id(id);  deallocate(x_solve)
    end do

    call init_id(id)
    id%ICNTL(21) = 1
    call attach_dense(id, n, nz, irn, jcn, A_trip, b)
    id%LSOL_loc = n
    allocate(id%SOL_loc(n))
    allocate(id%ISOL_loc(n))
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(21)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    ! Distributed solution (ICNTL(21)=1): each rank receives its own
    ! INFO(23) components in SOL_loc, with the global row indices in
    ! ISOL_loc; the union over ranks is the full solution. INFO(23) is
    ! the documented local-slice size (the struct's NSOL_loc field is NOT
    ! written by the solve). At np=1 the single rank owns all n
    ! components. max_rel_err_vec_z is max|Δ|/max|x_true| with a global
    ! denominator (x_true is full on every rank), so only the numerator —
    ! the max over ranks of each rank's local max modulus difference —
    ! needs a scalar MPI_MAX, reduced in double precision (see the
    ! d-variant for the rationale).
    block
        real(ep) :: local_num, denom
        real(8)  :: loc8, glob8
        integer  :: idx, gi
        local_num = 0.0_ep
        do idx = 1, id%INFO(23)
            gi = id%ISOL_loc(idx)
            local_num = max(local_num, &
                            abs(t2q_c(id%SOL_loc(idx)) - x_true(gi)))
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

    ! Centralized sparse RHS is a HOST-ONLY input: the solve driver
    ! rejects with INFO(1)=-22 any non-master rank whose RHS /
    ! RHS_SPARSE / IRHS_SPARSE / IRHS_PTR are associated. Only the host
    ! allocates them and the dense solution vector; slaves leave them
    ! null. The matrix (IRN/JCN/A) stays centralized on every rank.
    call init_id(id)
    id%ICNTL(20) = 1
    id%N    = n
    id%NNZ  = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = q2t_c(A_trip)
    id%NRHS    = 1
    if (myid == 0) then
        id%LRHS    = n
        id%NZ_RHS  = n
        allocate(id%RHS_SPARSE(n));   id%RHS_SPARSE  = q2t_c(b)
        allocate(id%IRHS_SPARSE(n))
        allocate(id%IRHS_PTR(2))
        do i = 1, n;  id%IRHS_SPARSE(i) = i;  end do
        id%IRHS_PTR(1) = 1
        id%IRHS_PTR(2) = n + 1
        allocate(id%RHS(n));  id%RHS = q2t_c(cmplx(0.0_ep, 0.0_ep, kind=ep))
    end if
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'ICNTL(20)=1 failed: INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    if (myid == 0) then
        allocate(x_solve(n));  x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
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
        type(zmumps_struc), intent(inout) :: id
        id%COMM = MPI_COMM_WORLD;  id%PAR = 1;  id%SYM = 0;  id%JOB = -1
        call target_xmumps(id)
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0
    end subroutine init_id

    subroutine attach_dense(id, n, nz, irn, jcn, A_trip, b)
        type(zmumps_struc), intent(inout) :: id
        integer,            intent(in)    :: n, nz, irn(:), jcn(:)
        complex(ep),        intent(in)    :: A_trip(:), b(:)
        id%N   = n
        id%NNZ = int(nz, kind=8)
        allocate(id%IRN(nz));  id%IRN = irn
        allocate(id%JCN(nz));  id%JCN = jcn
        allocate(id%A(nz));    id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));   id%RHS = q2t_c(b)
    end subroutine attach_dense

    subroutine end_id(id)
        type(zmumps_struc), intent(inout) :: id
        if (associated(id%IRN)) deallocate(id%IRN)
        if (associated(id%JCN)) deallocate(id%JCN)
        if (associated(id%A))   deallocate(id%A)
        if (associated(id%RHS)) deallocate(id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)
    end subroutine end_id

end program test_zmumps_icntl_io
