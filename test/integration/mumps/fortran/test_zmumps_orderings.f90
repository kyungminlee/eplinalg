! ICNTL(7) ordering coverage for ZMUMPS — mirror of test_dmumps_orderings.

program test_zmumps_orderings
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_sparse_hpd_problem_z
    use target_mumps,          only: target_name, zmumps_struc
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_solve6, mumps_check_solution, &
                                     mumps_end
    use mpi
    implicit none

    integer, parameter :: nx = 8, ny = 8
    integer, parameter :: orderings(*) = [0, 2, 3, 4, 5, 6, 7]
    integer            :: n, ierr, i, nz, ord
    complex(ep), allocatable :: x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_orderings', target_name)

    call gen_sparse_hpd_problem_z(nx, ny, n, x_true, b, irn, jcn, A_trip, nz, &
                                  seed = 24001)

    do i = 1, size(orderings)
        ord = orderings(i)

        write(label, '(a,i0)') 'icntl7=', ord
        call mumps_begin(id, MPI_COMM_WORLD, 0, ' ('//trim(label)//')')
        id%ICNTL(7) = ord

        call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
        call mumps_solve6(id, ' ('//trim(label)//')')

        call mumps_check_solution(id, x_true, trim(label))
        call mumps_end(id)
    end do

    deallocate(x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_orderings
