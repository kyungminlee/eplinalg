! Multiple-RHS coverage for ZMUMPS — mirror of test_dmumps_nrhs.

program test_zmumps_nrhs
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z, max_rel_err_mat_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer            :: ierr, nz, j
    complex(ep), allocatable :: A(:,:), x_true(:), b(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    complex(ep), allocatable :: B_multi(:,:), X_multi_true(:,:)
    complex(ep), allocatable :: x_single(:)
    real(ep)                 :: err, tol
    character(len=48)        :: label

    call MPI_INIT(ierr)
    call report_init('test_zmumps_nrhs', target_name)

    call gen_dense_problem_z(n, A, x_true, b, seed = 23001)
    call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
    tol = 16.0_ep * real(n, ep)**3 * target_eps

    call mumps_solve(n=n, nz=nz, nrhs=1, irn=irn, jcn=jcn, &
                     A_trip=A_trip, B_in=reshape(b, [n, 1]), &
                     X_out_buf=B_multi)
    allocate(x_single(n))
    x_single = B_multi(:, 1)
    err = max_rel_err_vec_z(x_single, x_true)
    call report_case('nrhs=1', err, tol)
    deallocate(B_multi, x_single)

    allocate(X_multi_true(n, 3), B_multi(n, 3))
    X_multi_true(:, 1) = x_true
    X_multi_true(:, 2) = cmplx(2.0_ep, 0.0_ep, kind=ep) * x_true &
                       - cmplx(1.0_ep, 0.0_ep, kind=ep)
    X_multi_true(:, 3) = cmplx(-0.5_ep, 0.0_ep, kind=ep) * x_true &
                       + cmplx(0.25_ep, 0.0_ep, kind=ep)
    B_multi = matmul(A, X_multi_true)

    block
        complex(ep), allocatable :: X_multi(:,:)
        call mumps_solve(n=n, nz=nz, nrhs=3, irn=irn, jcn=jcn, &
                         A_trip=A_trip, B_in=B_multi, X_out_buf=X_multi)
        err = max_rel_err_mat_z(X_multi, X_multi_true)
        call report_case('nrhs=3', err, tol)

        do j = 1, 3
            block
                complex(ep), allocatable :: single_col(:,:), b_col(:,:)
                allocate(b_col(n, 1))
                b_col(:, 1) = matmul(A, X_multi_true(:, j))
                call mumps_solve(n=n, nz=nz, nrhs=1, irn=irn, jcn=jcn, &
                                 A_trip=A_trip, B_in=b_col, X_out_buf=single_col)
                err = max_rel_err_vec_z(X_multi(:, j), single_col(:, 1))
                write(label, '(a,i0)') 'multi-vs-single-col=', j
                ! ZMUMPS multi-vs-single per-column is NOT bit-identical
                ! (the complex BLAS path through the solve driver
                ! reorders ops differently for NRHS=1 vs NRHS>1; D side
                ! happens to match because the real path is more
                ! straight-line). Use the standard n*eps tolerance.
                call report_case(trim(label), err, &
                                 16.0_ep * real(n, ep) * target_eps)
                deallocate(single_col, b_col)
            end block
        end do
        deallocate(X_multi)
    end block

    deallocate(A, x_true, b, irn, jcn, A_trip, B_multi, X_multi_true)
    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()

contains

    subroutine mumps_solve(n, nz, nrhs, irn, jcn, A_trip, B_in, X_out_buf)
        integer,     intent(in)               :: n, nz, nrhs, irn(:), jcn(:)
        complex(ep), intent(in)               :: A_trip(:), B_in(:,:)
        complex(ep), allocatable, intent(out) :: X_out_buf(:,:)
        type(zmumps_struc) :: idl
        integer :: i

        idl%COMM = MPI_COMM_WORLD;  idl%PAR = 1;  idl%SYM = 0;  idl%JOB = -1
        call target_xmumps(idl)
        idl%ICNTL(1) = -1; idl%ICNTL(2) = -1; idl%ICNTL(3) = -1; idl%ICNTL(4) = 0

        idl%N    = n
        idl%NNZ  = int(nz, kind=8)
        idl%NRHS = nrhs
        idl%LRHS = n
        allocate(idl%IRN(nz));  idl%IRN = irn
        allocate(idl%JCN(nz));  idl%JCN = jcn
        allocate(idl%A(nz));    idl%A   = q2t_c(A_trip)
        allocate(idl%RHS(n * nrhs))
        do i = 1, nrhs
            idl%RHS((i - 1) * n + 1 : i * n) = q2t_c(B_in(:, i))
        end do

        idl%JOB = 6
        call target_xmumps(idl)
        if (idl%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0)') 'NRHS=', nrhs, ' solve failed, INFOG(1)=', &
                idl%INFOG(1)
            error stop 1
        end if

        allocate(X_out_buf(n, nrhs))
        do i = 1, nrhs
            X_out_buf(:, i) = t2q_c(idl%RHS((i - 1) * n + 1 : i * n))
        end do

        deallocate(idl%IRN, idl%JCN, idl%A, idl%RHS)
        nullify(idl%IRN, idl%JCN, idl%A, idl%RHS)
        idl%JOB = -2
        call target_xmumps(idl)
    end subroutine mumps_solve

end program test_zmumps_nrhs
