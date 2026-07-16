! Complex mirror of test_dmumps_errors. See that file's header for
! rationale (D1 in tests/mumps/TODO.md).

program test_zmumps_errors
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, &
                                     report_finalize, report_check_status
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, zmumps_struc, target_xmumps, &
                                     q2t_c, &
                                     check_zmumps_input, &
                                     MIC_OK, MIC_BAD_N, MIC_BAD_NNZ, &
                                     MIC_BAD_IRN, MIC_BAD_JCN, MIC_SIZE_MISMATCH
    use mpi
    implicit none

    integer            :: ierr, nz, code, saved_n
    integer(8)         :: saved_nnz
    integer            :: saved_irn1, saved_jcn1
    integer, parameter :: n = 4
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc) :: id
    integer, pointer :: bad_irn(:)

    call MPI_INIT(ierr)
    call report_init('test_zmumps_errors', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 7919)
    call dense_to_triplet_z(A, irn, jcn, A_trip, nz)

    id%COMM = MPI_COMM_WORLD
    id%PAR  = 1
    id%SYM  = 0
    id%JOB  = -1
    call target_xmumps(id)
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
    allocate(id%A(nz));    id%A   = q2t_c(A_trip)
    allocate(id%RHS(n));   id%RHS = q2t_c(b)

    saved_n    = id%N
    saved_nnz  = id%NNZ
    saved_irn1 = id%IRN(1)
    saved_jcn1 = id%JCN(1)

    id%N = -1
    call check_zmumps_input(id, code)
    call expect('bad_n_neg',  code, MIC_BAD_N)
    id%N = 0
    call check_zmumps_input(id, code)
    call expect('bad_n_zero', code, MIC_BAD_N)
    id%N = saved_n

    id%NNZ = -3_8
    call check_zmumps_input(id, code)
    call expect('bad_nnz_neg', code, MIC_BAD_NNZ)
    id%NNZ = saved_nnz

    id%IRN(1) = saved_n + 5
    call check_zmumps_input(id, code)
    call expect('bad_irn_high', code, MIC_BAD_IRN)
    id%IRN(1) = 0
    call check_zmumps_input(id, code)
    call expect('bad_irn_zero', code, MIC_BAD_IRN)
    id%IRN(1) = saved_irn1

    id%JCN(1) = saved_n + 3
    call check_zmumps_input(id, code)
    call expect('bad_jcn_high', code, MIC_BAD_JCN)
    id%JCN(1) = 0
    call check_zmumps_input(id, code)
    call expect('bad_jcn_zero', code, MIC_BAD_JCN)
    id%JCN(1) = saved_jcn1

    bad_irn => id%IRN
    nullify(id%IRN)
    allocate(id%IRN(nz - 1));  id%IRN = bad_irn(1:nz-1)
    call check_zmumps_input(id, code)
    call expect('size_mismatch_irn', code, MIC_SIZE_MISMATCH)
    deallocate(id%IRN)
    id%IRN => bad_irn
    nullify(bad_irn)

    call check_zmumps_input(id, code)
    call expect('valid_ok', code, MIC_OK)
    id%JOB = 6
    call target_xmumps(id)
    if (id%INFOG(1) < 0) then
        write(*, '(a,i0)') 'JOB=6 on valid input failed, INFOG(1)=', id%INFOG(1)
        error stop 1
    end if

    deallocate(id%IRN, id%JCN, id%A, id%RHS)
    nullify(id%IRN, id%JCN, id%A, id%RHS)
    id%JOB = -2
    call target_xmumps(id)

    deallocate(A, x_true, b, irn, jcn, A_trip)

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine expect(label, got, want)
        character(len=*), intent(in) :: label
        integer, intent(in) :: got, want
        real(ep) :: err, tol
        if (got == want) then
            err = 0.0_ep
        else
            err = 1.0_ep
        end if
        tol = 0.5_ep
        call report_case(label, err, tol)
        if (got /= want) then
            write(*, '(a,a,a,i0,a,i0)') 'expect(', label, '): got ', got, &
                ' want ', want
        end if
    end subroutine expect

end program test_zmumps_errors
