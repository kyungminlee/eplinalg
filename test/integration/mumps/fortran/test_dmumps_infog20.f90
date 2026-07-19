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
    use target_mumps,          only: target_name, dmumps_struc
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_solve6, mumps_end
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

    call mumps_begin(id, MPI_COMM_WORLD, 0)
    call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
    call mumps_solve6(id)

    real_words = id%INFOG(20)

    block
        integer, parameter :: infog20_baseline = 1024   ! n=32 ⇒ n²
        write(label, '(a,i0,a,i0)') 'infog20=', real_words, ' n=', n
        err = abs(real(real_words - infog20_baseline, ep)) / &
              real(infog20_baseline, ep)
        tol = 0.05_ep
        call report_case(trim(label), err, tol)
    end block

    call mumps_end(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

end program test_dmumps_infog20
