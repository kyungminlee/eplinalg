! Genuine double-precision MUMPS Fortran-API solve test.
!
! Companion to c/genuine/test_genuine_dz_solve.c. That test drives the
! genuine solvers through the C bridges (dmumps_c / zmumps_c); this one
! calls the Fortran entry points DMUMPS / ZMUMPS directly on the
! DMUMPS_STRUC / ZMUMPS_STRUC derived types — the canonical upstream
! usage (mirrors external/MUMPS_5.8.2/examples/dsimpletest.F). It is the
! only per-stage ctest that exercises the genuine libdmumps / libzmumps
! Fortran API path, which the C bridge only reaches through the
! dmumps_f77_ shim; the derived-type layout, module linkage and Fortran
! calling convention are covered here and nowhere else.
!
! Unlike the migrated d/z-named Fortran tests under fortran/ (remapped
! onto q/x, e/y or m/w at extended precision), this calls the PRISTINE
! double-precision solver, so residuals come back at plain-double
! accuracy (~1e-16), not the migrated width.
!
! Centralized-input MUMPS: matrix + RHS supplied on the host only, the
! centralized solution returned in id%RHS on the host only. DMUMPS /
! ZMUMPS calls stay collective on every rank; population and the
! residual check are host-guarded, so the test is correct under
! mpiexec -n>1 as well as MPISEQ. The Fortran communicator is passed
! directly via id%COMM (no MPI_Comm_c2f mapping, unlike the C bridge).

program test_genuine_dz_solve_f
    use mpi
    implicit none
    integer :: myid, ierr
    logical :: is_host, ok_r, ok_z

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    is_host = (myid == 0)

    call solve_real(is_host, ok_r)
    call solve_complex(is_host, ok_z)

    call MPI_FINALIZE(ierr)
    if (is_host) then
        if (ok_r .and. ok_z) then
            stop 0
        else
            error stop 1
        end if
    end if

contains

    ! ── real solve via genuine DMUMPS ────────────────────────────────
    subroutine solve_real(is_host, ok)
        logical, intent(in)  :: is_host
        logical, intent(out) :: ok
        include 'dmumps_struc.h'
        type(dmumps_struc) :: id
        integer, parameter  :: n = 4
        double precision :: Ar(n,n), xt(n), mr, den, tol
        integer :: i, j, k

        ! Diagonally-dominant 4x4 (column-major reshape of the same
        ! matrix c/genuine/test_genuine_dz_solve.c uses).
        Ar = reshape([ &
             5.0d0,  0.3d0, -0.4d0,  0.1d0, &
             0.5d0, -6.0d0,  0.2d0, -0.3d0, &
            -0.25d0, 0.5d0,  7.0d0,  0.4d0, &
             0.1d0,  0.2d0, -0.3d0, -8.0d0 ], [n, n])
        xt = [1.0d0, -2.0d0, 3.0d0, -4.0d0]

        id%COMM = MPI_COMM_WORLD
        id%PAR = 1; id%SYM = 0; id%JOB = -1
        call DMUMPS(id)
        if (id%INFOG(1) < 0) then; ok = .false.; return; end if
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0

        if (is_host) then
            id%N   = n
            id%NNZ = int(n*n, kind=8)
            allocate(id%IRN(n*n), id%JCN(n*n), id%A(n*n), id%RHS(n))
            k = 0
            do j = 1, n
                do i = 1, n
                    k = k + 1
                    id%IRN(k) = i; id%JCN(k) = j; id%A(k) = Ar(i,j)
                end do
            end do
            do i = 1, n
                id%RHS(i) = sum(Ar(i,:) * xt)
            end do
        end if

        id%JOB = 6
        call DMUMPS(id)
        if (id%INFOG(1) < 0) then; ok = .false.; return; end if

        ok = .true.
        if (is_host) then
            den = maxval(abs(xt))
            mr  = maxval(abs(id%RHS - xt))
            if (den > 0.0d0) mr = mr / den
            tol = 16.0d0 * real(n, 8)**3 * epsilon(1.0d0)
            ok  = (mr <= tol)
            write(*, '(a,es13.6,a)') &
                '  test_genuine_dz_solve_f [real-n=4]    max_rel_err=', mr, &
                merge('  PASS', '  FAIL', ok)
            deallocate(id%IRN, id%JCN, id%A, id%RHS)
        end if

        id%JOB = -2
        call DMUMPS(id)
    end subroutine solve_real

    ! ── complex solve via genuine ZMUMPS ─────────────────────────────
    subroutine solve_complex(is_host, ok)
        logical, intent(in)  :: is_host
        logical, intent(out) :: ok
        include 'zmumps_struc.h'
        type(zmumps_struc) :: id
        integer, parameter  :: n = 4
        double precision :: Ar(n,n), Ai(n,n), mr, den, tol
        complex(kind=8)  :: A(n,n), xt(n)
        integer :: i, j, k

        Ar = reshape([ &
             5.0d0,  0.3d0, -0.4d0,  0.1d0, &
             0.5d0, -6.0d0,  0.2d0, -0.3d0, &
            -0.25d0, 0.5d0,  7.0d0,  0.4d0, &
             0.1d0,  0.2d0, -0.3d0, -8.0d0 ], [n, n])
        Ai = reshape([ &
             0.2d0,  0.03d0, 0.0d0,  0.01d0, &
             0.05d0, 0.3d0,  0.02d0, 0.0d0,  &
             0.0d0, -0.05d0,-0.4d0,  0.04d0, &
             0.01d0, 0.0d0,  0.03d0, 0.5d0 ], [n, n])
        A  = cmplx(Ar, Ai, kind=8)
        xt = [ cmplx( 1.0d0,  1.0d0, kind=8), cmplx(-2.0d0, 0.5d0, kind=8), &
               cmplx( 3.0d0, -1.0d0, kind=8), cmplx(-4.0d0, 2.0d0, kind=8) ]

        id%COMM = MPI_COMM_WORLD
        id%PAR = 1; id%SYM = 0; id%JOB = -1
        call ZMUMPS(id)
        if (id%INFOG(1) < 0) then; ok = .false.; return; end if
        id%ICNTL(1) = -1; id%ICNTL(2) = -1; id%ICNTL(3) = -1; id%ICNTL(4) = 0

        if (is_host) then
            id%N   = n
            id%NNZ = int(n*n, kind=8)
            allocate(id%IRN(n*n), id%JCN(n*n), id%A(n*n), id%RHS(n))
            k = 0
            do j = 1, n
                do i = 1, n
                    k = k + 1
                    id%IRN(k) = i; id%JCN(k) = j; id%A(k) = A(i,j)
                end do
            end do
            do i = 1, n
                id%RHS(i) = sum(A(i,:) * xt)
            end do
        end if

        id%JOB = 6
        call ZMUMPS(id)
        if (id%INFOG(1) < 0) then; ok = .false.; return; end if

        ok = .true.
        if (is_host) then
            den = maxval(abs(xt))
            mr  = maxval(abs(id%RHS - xt))
            if (den > 0.0d0) mr = mr / den
            tol = 16.0d0 * real(n, 8)**3 * epsilon(1.0d0)
            ok  = (mr <= tol)
            write(*, '(a,es13.6,a)') &
                '  test_genuine_dz_solve_f [complex-n=4] max_rel_err=', mr, &
                merge('  PASS', '  FAIL', ok)
            deallocate(id%IRN, id%JCN, id%A, id%RHS)
        end if

        id%JOB = -2
        call ZMUMPS(id)
    end subroutine solve_complex

end program test_genuine_dz_solve_f
