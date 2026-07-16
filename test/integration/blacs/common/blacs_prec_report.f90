! Module name is blacs_prec_report (not prec_report) to avoid colliding
! with tests/blas/common/prec_report and pbblas_prec_report — every
! framework writes its .mod files into the same ${PROJECT_BINARY_DIR}/fmod
! directory, so the module name must be unique.
!
! BLACS tests are exact (no roundoff) — callers pass tolerance=0.0_ep
! and supply err=0.0_ep on success or err=huge(0.0_ep) on mismatch.
! report_init/finalize are collective; only rank 0 opens/closes the
! JSON file but the failure flag is broadcast so every rank exits with
! the same status (so ctest sees a uniform pass/fail).
module blacs_prec_report
    use prec_kinds, only: ep
    use mpi
    implicit none
    private
    public :: report_init, report_case, report_finalize

    integer :: unit_save = -1
    integer :: case_count = 0
    integer :: rank_save = -1
    character(len=:), allocatable :: routine_save
    logical :: any_failure = .false.

contains

    subroutine report_init(routine, target_name, rank)
        character(len=*), intent(in) :: routine, target_name
        integer,          intent(in) :: rank
        character(len=:), allocatable :: filename
        integer :: ios

        routine_save = trim(routine)
        rank_save = rank
        any_failure = .false.
        case_count = 0

        if (rank /= 0) return

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

        if (rank_save /= 0) return

        if (max_rel > 0.0_ep) then
            digits = -log10(max_rel)
        else
            digits = 99.0_ep
        end if
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
            '] err=', max_rel, '  digits=', digits, &
            merge(' PASS', ' FAIL', passed)
    end subroutine report_case

    subroutine report_finalize()
        integer :: failure_flag, ierr

        if (rank_save == 0) then
            write(unit_save, '(a)') '  ]'
            write(unit_save, '(a)') '}'
            close(unit_save)
            unit_save = -1
        end if

        failure_flag = 0
        if (any_failure) failure_flag = 1
        call mpi_bcast(failure_flag, 1, mpi_integer, 0, mpi_comm_world, ierr)
        if (failure_flag /= 0) error stop 1
    end subroutine report_finalize

end module blacs_prec_report
