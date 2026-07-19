! Shared JSON precision-report module for the integration-test suites.
!
! One preprocessed source replaces the per-suite copies that used to
! live in tests/{blas,xblas,lapack,ptzblas,pblas,pbblas,scalapack,
! blacs,mumps}/common. Each suite compiles this file into its own
! Fortran_MODULE_DIRECTORY, selecting a variant with cpp macros:
!
!   (no macros)               serial reporter (ptzblas).
!   PREC_REPORT_STRICT_EXACT  serial reporter with the bit-exact
!                             kind16 gate (blas, xblas, lapack).
!   PREC_REPORT_MPI           MPI reporter: rank 0 writes the JSON and
!                             report_finalize broadcasts the pass/fail
!                             flag so every rank exits with the same
!                             status (pblas, pbblas, scalapack, blacs).
!   PREC_REPORT_MPI + PREC_REPORT_DEFERRED_STATUS
!                             MPI reporter whose report_finalize only
!                             closes the JSON; drivers run MPI_FINALIZE
!                             and then report_check_status() to gate
!                             the exit code (mumps — including the
!                             libmpiseq _seq variant, which resolves
!                             MPI_COMM_RANK against the stub library;
!                             no collective is referenced in this
!                             configuration).
!
! report_init's `rank` argument is optional: suites that track the
! rank themselves pass it (the pblas-family call shape), the others
! omit it and the module queries MPI_COMM_RANK itself (MPI variants)
! or reports from rank 0 unconditionally (serial variants).
!
! The historical per-suite module renames (pblas_prec_report and
! friends) survive as one-line re-export shims next to this file so
! existing test drivers keep their `use` lines unchanged.
module prec_report
    use prec_kinds, only: ep
#ifdef PREC_REPORT_MPI
    use mpi
#endif
    implicit none
    private
    public :: report_init, report_case, report_finalize, report_check_status

    integer :: unit_save = -1
    integer :: case_count = 0
    integer :: rank_save = -1
    character(len=:), allocatable :: routine_save
    logical :: any_failure = .false.
#ifdef PREC_REPORT_STRICT_EXACT
    ! Bit-exact mode for kind16: the migrated routine and the
    ! quad-promoted Netlib reference compute the same serial algorithm
    ! at REAL(KIND=16) — any non-zero divergence is a real migration
    ! bug, even one that lands inside the per-case rounding budget.
    ! (Distributed and sparse suites never enable this: their
    ! reductions and factorizations reorder the rounding.)
    logical :: strict_exact = .false.
#endif

contains

    ! report_init is called by every rank (collective in the MPI
    ! variants). Only rank 0 opens the JSON file; the other ranks just
    ! record their rank so report_case / report_finalize stay silent
    ! off-host.
    subroutine report_init(routine, target_name, rank)
        character(len=*), intent(in) :: routine, target_name
        integer, intent(in), optional :: rank
        character(len=:), allocatable :: filename
        integer :: ios
#ifdef PREC_REPORT_MPI
        integer :: ierr, myid
#endif

        if (present(rank)) then
            rank_save = rank
        else
#ifdef PREC_REPORT_MPI
            ! Every MPI driver calls MPI_INIT before report_init, so
            ! the communicator is live here. (Under the libmpiseq _seq
            ! variant MPI_COMM_RANK returns 0, so rank 0 still hosts.)
            call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
            rank_save = myid
#else
            rank_save = 0
#endif
        end if

        routine_save = trim(routine)
        any_failure = .false.
        case_count = 0
#ifdef PREC_REPORT_STRICT_EXACT
        strict_exact = (trim(target_name) == 'kind16')
#endif

        if (rank_save /= 0) return

        filename = trim(routine) // '.' // trim(target_name) // '.json'
        open(newunit=unit_save, file=filename, status='replace', &
             action='write', iostat=ios)
        if (ios /= 0) then
            write(*, '(a,a)') 'report_init: cannot open ', filename
            error stop 1
        end if

        write(unit_save, '(a)')      '{'
        write(unit_save, '(a,a,a)')  '  "routine": "', trim(routine), '",'
        write(unit_save, '(a,a,a)')  '  "target":  "', trim(target_name), '",'
        write(unit_save, '(a)')      '  "cases": ['
    end subroutine report_init

    ! report_case runs on rank 0 only — drivers that guard their
    ! comparison block on their own rank never enter here off-host;
    ! the rank_save check is defensive.
    subroutine report_case(case_label, max_rel, tol)
        character(len=*), intent(in) :: case_label
        real(ep), intent(in) :: max_rel, tol
        character(len=32) :: relbuf, tolbuf, digbuf
        character(len=5)  :: passbuf
        logical  :: passed
        real(ep) :: digits

        if (rank_save /= 0) return

        if (max_rel > 0.0_ep) then
            digits = -log10(max_rel)
        else
            digits = 99.0_ep
        end if
#ifdef PREC_REPORT_STRICT_EXACT
        if (strict_exact) then
            ! kind16 strict floor: ≥ 28 decimal digits of agreement
            ! against the quad reference. Quad arithmetic gives ~33
            ! significant digits, so 28 leaves ~5-ULP headroom for
            ! intrinsic-call differences (sqrt, divisions) and the
            ! ULPs accumulated by self-residual computations in the
            ! eigenvalue/SVD drivers — without ever accepting the
            ! kind of 1–5 % bug a precision-cast slip-up produces.
            passed = (max_rel == 0.0_ep) .or. (digits >= 28.0_ep)
        else
            passed = (max_rel <= tol)
        end if
#else
        passed = (max_rel <= tol)
#endif
        if (.not. passed) any_failure = .true.

        write(relbuf, '(es15.6e3)') max_rel
        write(tolbuf, '(es15.6e3)') tol
        write(digbuf, '(f6.2)')     digits
        passbuf = merge('true ', 'false', passed)

        if (case_count > 0) write(unit_save, '(a)') '    ,'
        case_count = case_count + 1

        write(unit_save, '(a)')     '    {'
        write(unit_save, '(a,a,a)') '      "case":        "', trim(case_label), '",'
        write(unit_save, '(a,a,a)') '      "max_rel_err": ', trim(adjustl(relbuf)), ','
        write(unit_save, '(a,a,a)') '      "tolerance":   ', trim(adjustl(tolbuf)), ','
        write(unit_save, '(a,a,a)') '      "digits":      ', trim(adjustl(digbuf)), ','
        write(unit_save, '(a,a)')   '      "passed":      ', trim(passbuf)
        write(unit_save, '(a)')     '    }'
        ! Flush after every case so a downstream crash leaves an
        ! inspectable JSON tail rather than buffered nothingness.
        flush(unit_save)

        ! Echo to stdout for ctest --output-on-failure / human reading.
        write(*, '(a,a,a,a,a,es12.4,a,f6.2,a)') &
            '  ', trim(routine_save), ' [', trim(case_label), &
            '] max_rel_err=', max_rel, '  digits=', digits, &
            merge(' PASS', ' FAIL', passed)
    end subroutine report_case

    subroutine report_finalize()
#if defined(PREC_REPORT_MPI) && !defined(PREC_REPORT_DEFERRED_STATUS)
        integer :: failure_flag, ierr
#endif

        if (rank_save == 0) then
            write(unit_save, '(a)') '  ]'
            write(unit_save, '(a)') '}'
            flush(unit_save)
            close(unit_save)
            unit_save = -1
        end if

#ifdef PREC_REPORT_DEFERRED_STATUS
        ! Deferred status: no error stop here even when cases failed —
        ! drivers must call MPI_FINALIZE first and then invoke
        ! report_check_status() to halt the program. Splitting the two
        ! avoids skipping MPI cleanup on failure (an `error stop 1`
        ! here would leave MPI in a state the runtime prints warnings
        ! about).
#elif defined(PREC_REPORT_MPI)
        ! report_finalize is collective — every rank calls it. Rank 0
        ! broadcasts the pass/fail flag so every rank exits with the
        ! same status (needed so ctest sees a uniform failure/success
        ! across the whole mpirun invocation).
        failure_flag = 0
        if (any_failure) failure_flag = 1
        call mpi_bcast(failure_flag, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (failure_flag /= 0) error stop 1
#else
        ! Non-zero exit code if any case failed — CTest sees pass/fail.
        if (any_failure) error stop 1
#endif
    end subroutine report_finalize

    ! Halts with a non-zero exit code when any case failed. The
    ! deferred-status drivers call this after MPI_FINALIZE; for the
    ! other variants report_finalize has already gated the exit and
    ! this is a no-op safety net.
    subroutine report_check_status()
        if (rank_save == 0 .and. any_failure) error stop 1
    end subroutine report_check_status

end module prec_report
