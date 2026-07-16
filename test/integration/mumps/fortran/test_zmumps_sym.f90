! Symmetry-flag coverage for ZMUMPS:
!   SYM = 0 — general unsymmetric complex (LU)
!   SYM = 1 — complex symmetric positive-definite. ZMUMPS implements
!             SYM=1 as A = L*L^T (note: NOT L*L^H — Cholesky-style on
!             a complex-SYMMETRIC matrix, not Hermitian). Easiest
!             construction: a real SPD matrix cast to complex with
!             zero imaginary parts, which trivially satisfies A=A^T
!             and has positive eigenvalues.
!   SYM = 2 — complex symmetric (NOT Hermitian — A = A^T, not A^H)

program test_zmumps_sym
    use prec_kinds,            only: ep
    use prec_report,           only: report_init, report_case, report_finalize, report_check_status
    use compare,               only: max_rel_err_vec_z
    use test_data_mumps,       only: gen_dense_problem_z, dense_to_triplet_z, &
                                     dense_to_sym_triplet_z, gen_spd_dense_problem
    use target_mumps,          only: target_name, target_eps, &
                                     zmumps_struc, target_xmumps, &
                                     q2t_c, t2q_c
    use mpi
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: syms(*) = [0, 1, 2]
    integer            :: ierr, k, nz, sym
    complex(ep), allocatable :: A(:,:), A_sym(:,:), x_true(:), b(:), x_solve(:)
    real(ep),    allocatable :: A_real(:,:), x_true_real(:), b_real(:)
    integer,     allocatable :: irn(:), jcn(:)
    complex(ep), allocatable :: A_trip(:)
    type(zmumps_struc)       :: id
    real(ep)                 :: err, tol
    character(len=48)        :: label
    integer :: i, j

    call MPI_INIT(ierr)
    call report_init('test_zmumps_sym', target_name)

    do k = 1, size(syms)
        sym = syms(k)

        if (sym == 0) then
            call gen_dense_problem_z(n, A, x_true, b, seed = 8001 + k)
            call dense_to_triplet_z (A, irn, jcn, A_trip, nz)
        else if (sym == 1) then
            ! Complex SPD = real SPD cast to complex. A = X*X^T + n*I
            ! is real symmetric positive-definite; embedding it in a
            ! complex matrix with zero imaginary part keeps A = A^T
            ! AND positive-definite for ZMUMPS's L*L^T factorization.
            call gen_spd_dense_problem(n, A_real, x_true_real, b_real, &
                                       seed = 8021 + k)
            allocate(A(n, n), x_true(n), b(n))
            A      = cmplx(A_real,      0.0_ep, kind=ep)
            x_true = cmplx(x_true_real, 0.0_ep, kind=ep)
            b      = cmplx(b_real,      0.0_ep, kind=ep)
            deallocate(A_real, x_true_real, b_real)
            call dense_to_sym_triplet_z(A, irn, jcn, A_trip, nz)
        else
            ! SYM=2: complex-symmetric A = (R + R^T)/2 with diagonal
            ! boost. Indefinite is fine; only structure A=A^T matters.
            call gen_dense_problem_z(n, A, x_true, b, seed = 8011 + k)
            allocate(A_sym(n, n))
            A_sym = 0.5_ep * (A + transpose(A))
            do i = 1, n
                A_sym(i, i) = A_sym(i, i) + cmplx(real(n, ep), 0.0_ep, kind=ep)
            end do
            ! recompute b for the symmetrized A
            b = matmul(A_sym, x_true)
            call move_alloc(A_sym, A)
            call dense_to_sym_triplet_z(A, irn, jcn, A_trip, nz)
        end if

        id%COMM = MPI_COMM_WORLD
        id%PAR  = 1
        id%SYM  = sym
        id%JOB  = -1
        call target_xmumps(id)

        id%ICNTL(1) = -1
        id%ICNTL(2) = -1
        id%ICNTL(3) = -1
        id%ICNTL(4) = 0

        id%N    = n
        id%NNZ  = int(nz, kind=8)
        allocate(id%IRN(nz));    id%IRN = irn
        allocate(id%JCN(nz));    id%JCN = jcn
        allocate(id%A(nz));      id%A   = q2t_c(A_trip)
        allocate(id%RHS(n));     id%RHS = q2t_c(b)

        id%JOB = 6
        call target_xmumps(id)
        if (id%INFOG(1) < 0) then
            write(*, '(a,i0,a,i0,a,i0)') 'JOB=6 (sym=', sym, &
                ') failed, INFOG(1)=', id%INFOG(1), &
                ', INFOG(2)=', id%INFOG(2)
            error stop 1
        end if

        allocate(x_solve(n))
        x_solve = t2q_c(id%RHS)
        err = max_rel_err_vec_z(x_solve, x_true)
        tol = 16.0_ep * real(n, ep)**3 * target_eps
        write(label, '(a,i0)') 'sym=', sym
        call report_case(trim(label), err, tol)

        deallocate(id%IRN, id%JCN, id%A, id%RHS)
        nullify(id%IRN, id%JCN, id%A, id%RHS)
        id%JOB = -2
        call target_xmumps(id)

        deallocate(A, x_true, b, x_solve, irn, jcn, A_trip)
    end do

    call report_finalize()
    call MPI_FINALIZE(ierr)
    call report_check_status()
end program test_zmumps_sym
