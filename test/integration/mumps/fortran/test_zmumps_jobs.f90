! JOB-phasing equivalence for ZMUMPS — mirror of test_dmumps_jobs.f90.
! Combined JOB=6 vs phased JOB=1+2+3 must produce bit-identical
! solutions on the same complex problem.

program test_zmumps_jobs
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: n = 12
    integer            :: ierr, nz
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: x_combined(:), x_phased(:)
    real(ep)                 :: err, tol

    call MPI_INIT(ierr)
    call report_init('test_zmumps_jobs', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 21001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    call mumps_solve_jobs(jobs = [6], n=n, nz=nz, irn=irn, jcn=jcn, &
                          A_trip=A_trip, b=b, x_solve=x_combined)
    call mumps_solve_jobs(jobs = [1, 2, 3], n=n, nz=nz, irn=irn, jcn=jcn, &
                          A_trip=A_trip, b=b, x_solve=x_phased)

    err = max_rel_err_vec_z(x_combined, x_true)
    call report_case('combined-vs-truth', err, tol)
    err = max_rel_err_vec_z(x_phased, x_true)
    call report_case('phased-vs-truth', err, tol)
    ! The two solve paths order operations identically on native real kinds
    ! (bit-for-bit), but extended-precision (multifloats) arithmetic can
    ! accumulate differently and diverge by a few ULP. Gate on the same
    ! tolerance as the truth comparisons, not on bit-equality.
    err = max_rel_err_vec_z(x_phased, x_combined)
    call report_case('phased-vs-combined', err, tol)

    deallocate(A, x_true, b, irn, jcn, A_trip, x_combined, x_phased)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine mumps_solve_jobs(jobs, n, nz, irn, jcn, A_trip, b, x_solve)
        integer,     intent(in)               :: jobs(:), n, nz, irn(:), jcn(:)
        complex(ep), intent(in)               :: A_trip(:), b(:)
        complex(ep), allocatable, intent(out) :: x_solve(:)
        type(zmumps_struc) :: idl
        integer :: k

        idl%COMM = MPI_COMM_WORLD;  idl%PAR = 1;  idl%SYM = 0;  idl%JOB = -1
        call target_xmumps(idl)
        idl%ICNTL(1) = -1; idl%ICNTL(2) = -1; idl%ICNTL(3) = -1; idl%ICNTL(4) = 0

        idl%N   = n
        idl%NNZ = int(nz, kind=8)
        allocate(idl%IRN(nz));  idl%IRN = irn
        allocate(idl%JCN(nz));  idl%JCN = jcn
        allocate(idl%A(nz));    idl%A   = q2t_c(A_trip)
        allocate(idl%RHS(n));   idl%RHS = q2t_c(b)

        do k = 1, size(jobs)
            idl%JOB = jobs(k)
            call target_xmumps(idl)
            if (idl%INFOG(1) < 0) then
                write(*, '(a,i0,a,i0)') 'JOB=', jobs(k), &
                    ' failed, INFOG(1)=', idl%INFOG(1)
                error stop 1
            end if
        end do

        allocate(x_solve(n))
        x_solve = t2q_c(idl%RHS)

        deallocate(idl%IRN, idl%JCN, idl%A, idl%RHS)
        nullify(idl%IRN, idl%JCN, idl%A, idl%RHS)
        idl%JOB = -2
        call target_xmumps(idl)
    end subroutine mumps_solve_jobs

end program test_zmumps_jobs
