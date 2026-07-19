! Basic smoke + tolerance baseline for the migrated ZMUMPS via the
! kind16 ${LIB_PREFIX}mumps archive (called as `xmumps` here through
! the target_mumps wrapper). Mirrors test_dmumps_basic.f90 but with
! complex-valued coefficient and right-hand side.

program test_zmumps_basic
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, zmumps_struc
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_solve6, mumps_check_solution, &
                                     mumps_end
    use mpi
    implicit none

    integer, parameter :: ns(*) = [8, 32]
    integer            :: ierr, i, n, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_basic', target_name)

    do i = 1, size(ns)
        n = ns(i)
        call gen_dense_problem_z(n, A, x_true, b, seed = 7001 + 31 * i)
        call dense_to_triplet_z (A, irn, jcn, A_trip, nz)

        call mumps_begin(id, MPI_COMM_WORLD, 0)
        call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
        call mumps_solve6(id)

        write(label, '(a,i0)') 'n=', n
        call mumps_check_solution(id, x_true, trim(label))
        call mumps_end(id)

        deallocate(A, x_true, b, irn, jcn, A_trip)
    end do

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_basic
