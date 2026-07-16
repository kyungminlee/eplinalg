! Cross-target sanity check on INFOG(20) — peak number of
! real-precision words used during factorization. The element
! count is structurally precision-independent for a fixed
! (n, SYM, JOB) problem (the multifrontal tree shape depends on
! matrix structure, not precision), so any drift across kind10 /
! kind16 / multifloats points at the hand-written byte-accounting
! overrides for `mumps_memory_mod` / `mumps_lr_stats`.
!
! Captured baseline (n=32, SYM=0, JOB=6, seed=4099): INFOG(20) = 1024
! on all three targets, bit-exact. A ±5% tolerance against that
! baseline catches override-sizing drift well below the gross
! mis-sizing threshold the original structural [n², 50·n²] window
! covered.

program test_dmumps_infog20
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, &
                                     report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem, dense_to_triplet
    use target_mumps,          only: target_name, dmumps_struc, target_qmumps, &
                                     q2t_r
    use mpi
    implicit none

    integer, parameter :: n = 32
    integer            :: ierr, nz
    real(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc) :: id
    integer            :: real_words
    real(ep)           :: err, tol
    character(len=64)  :: label

    call MPI_INIT(ierr)
    call report_init('test_dmumps_infog20', target_name)

    call gen_dense_problem(n, A, x_true, b, seed = 4099)
    call dense_to_triplet(A, irn, jcn, A_trip, nz)

    id%COMM = MPI_COMM_WORLD
    id%PAR  = 1
    id%SYM  = 0
    id%JOB  = -1
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'JOB=-1 failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if
    id%ICNTL(1:3) = -1
    id%ICNTL(4)   = 0

    id%N    = n
    id%NNZ  = int(nz, kind=8)
    allocate(id%IRN(nz));  id%IRN = irn
    allocate(id%JCN(nz));  id%JCN = jcn
    allocate(id%A(nz));    id%A   = q2t_r(A_trip)
    allocate(id%RHS(n));   id%RHS = q2t_r(b)

    id%JOB = 6
    call target_qmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0,a,i0)') 'JOB=6 failed, INFOG(1)=', &
            id%INFOG(1), ', INFOG(2)=', id%INFOG(2)
        error stop 1
    end if

    real_words = id%INFOG(20)

    block
        integer, parameter :: infog20_baseline = 1024   ! n=32 ⇒ n²
        write(label, '(a,i0,a,i0)') 'infog20=', real_words, ' n=', n
        err = abs(real(real_words - infog20_baseline, ep)) / &
              real(infog20_baseline, ep)
        tol = 0.05_ep
        call report_case(trim(label), err, tol)
    end block

    deallocate(id%IRN, id%JCN, id%A, id%RHS)
    nullify(id%IRN, id%JCN, id%A, id%RHS)
    id%JOB = -2
    call target_qmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

end program test_dmumps_infog20
