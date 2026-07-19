! ICNTL(7) ordering coverage for DMUMPS:
!   0 — AMD (approximate minimum degree)
!   2 — AMF (approximate minimum fill)
!   3 — Scotch (nested dissection; vendored + privately namespaced)
!   4 — PORD (nested dissection; ships in-tree with MUMPS)
!   5 — METIS (nested dissection; vendored + privately namespaced)
!   6 — QAMD (with extra quasi-dense rows)
!   7 — automatic (let MUMPS choose)
!
! Each ordering is a different graph-partition strategy fed into the
! analysis phase; they alter the elimination tree and thus the exact
! floating-point order of the factorization, but should produce the
! same numeric solution to within a fairly tight tolerance against
! the quad-precision ground truth.
!
! The problem is a SPARSE 2D-Laplacian SPD system (not the dense
! generators the other tests use): nested-dissection orderings such as
! PORD need a graph with real vertex separators — on a dense/complete
! graph they abort with "no valid number of stages in multisector".

program test_dmumps_orderings
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_finalize, report_check_status
    use test_data_mumps,       only: gen_sparse_spd_problem
    use target_mumps,          only: target_name, dmumps_struc
    use mumps_lifecycle,       only: mumps_begin, mumps_load_triplet, &
                                     mumps_solve6, mumps_check_solution, &
                                     mumps_end
    use mpi
    implicit none

    integer, parameter :: nx = 8, ny = 8
    integer, parameter :: orderings(*) = [0, 2, 3, 4, 5, 6, 7]
    integer            :: n, ierr, i, nz, ord
    real(ep), allocatable :: x_true(:), b(:)
    integer,  allocatable :: irn(:), jcn(:)
    real(ep), allocatable :: A_trip(:)
    type(dmumps_struc)    :: id
    character(len=48)     :: label

    call MPI_INIT(ierr)
    call report_init('test_dmumps_orderings', target_name)

    call gen_sparse_spd_problem(nx, ny, n, x_true, b, irn, jcn, A_trip, nz, &
                                seed = 5001)

    do i = 1, size(orderings)
        ord = orderings(i)

        write(label, '(a,i0)') 'icntl7=', ord
        call mumps_begin(id, MPI_COMM_WORLD, 0, ' ('//trim(label)//')')
        id%ICNTL(7) = ord    ! analysis ordering selector

        call mumps_load_triplet(id, n, nz, irn, jcn, A_trip, b)
        call mumps_solve6(id, ' ('//trim(label)//')')

        call mumps_check_solution(id, x_true, trim(label))
        call mumps_end(id)
    end do

    deallocate(x_true, b, irn, jcn, A_trip)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_dmumps_orderings
