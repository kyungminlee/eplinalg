module prec_report
    use prec_kinds, only: ep
    use mpi
    implicit none
    private
    public :: report_init, report_case, report_finalize, report_check_status

    integer :: unit_save = -1
    integer :: case_count = 0
    character(len=:), allocatable :: routine_save
    logical :: any_failure = .false.
    ! MUMPS centralizes the solution on the host (rank 0); non-host ranks
    ! never hold a checkable result. So all reporting — the JSON file I/O,
    ! the per-case check, and the final error stop — is confined to the
    ! host, decided once here from the MPI rank. This lets an otherwise
    ! rank-agnostic driver run correctly under ``mpiexec -n N``: the
    ! collective MUMPS calls stay on every rank while only rank 0 writes
    ! the report and gates the exit status. At np=1 the host is the only
    ! rank, so behaviour is identical to the single-rank path. Drivers
    ! that already wrap their report_* calls in ``if (is_host)`` are
    ! unaffected (they never enter these routines off-host).
    logical :: is_host = .true.

contains

    subroutine report_init(routine, target_name)
        character(len=*), intent(in) :: routine, target_name
        character(len=:), allocatable :: filename
        integer :: ios, ierr, myid

        ! Every MUMPS driver calls MPI_INIT before report_init, so the
        ! communicator is live here. (Under the libmpiseq _seq variant
        ! MPI_COMM_RANK returns 0, so is_host stays .true..)
        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
        is_host = (myid == 0)

        any_failure = .false.
        case_count = 0
        if (.not. is_host) return

        routine_save = trim(routine)
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

    subroutine report_case(case_label, max_rel, tol)
        character(len=*), intent(in) :: case_label
        real(ep), intent(in) :: max_rel, tol
        character(len=32) :: relbuf, tolbuf, digbuf
        character(len=5)  :: passbuf
        logical  :: passed
        real(ep) :: digits

        if (.not. is_host) return

        if (max_rel > 0.0_ep) then
            digits = -log10(max_rel)
        else
            digits = 99.0_ep
        end if
        ! No strict_exact mode for MUMPS — the reference (reflapack_quad
        ! qdgesv) is dense LU at REAL(KIND=16) while the migrated
        ! ${LIB_PREFIX}mumps does sparse direct factorization. Same
        ! mathematical result up to algorithmic rounding-order
        ! differences, so the per-case tolerance is the only valid gate.
        passed = (max_rel <= tol)
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

        write(*, '(a,a,a,a,a,es12.4,a,f6.2,a)') &
            '  ', trim(routine_save), ' [', trim(case_label), &
            '] max_rel_err=', max_rel, '  digits=', digits, &
            merge(' PASS', ' FAIL', passed)
    end subroutine report_case

    subroutine report_finalize()
        ! Closes the JSON report. Does NOT raise an error when cases
        ! failed — callers must call MPI_FINALIZE first then invoke
        ! report_check_status() to halt the program. Splitting the two
        ! avoids skipping MPI cleanup on failure (the previous combined
        ! routine called `error stop 1` here, leaving MPI in a state
        ! that the runtime would print warnings about).
        if (.not. is_host) return
        write(unit_save, '(a)') '  ]'
        write(unit_save, '(a)') '}'
        close(unit_save)
        unit_save = -1
    end subroutine report_finalize

    subroutine report_check_status()
        if (is_host .and. any_failure) error stop 1
    end subroutine report_check_status

end module prec_report
