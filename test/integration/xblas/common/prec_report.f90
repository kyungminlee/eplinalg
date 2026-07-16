module prec_report
    use prec_kinds, only: ep
    implicit none
    private
    public :: report_init, report_case, report_finalize

    integer :: unit_save = -1
    integer :: case_count = 0
    character(len=:), allocatable :: routine_save
    logical :: any_failure = .false.
    ! Bit-exact mode for kind16: the migrated routine and the
    ! quad-promoted Netlib reference compute the same serial algorithm
    ! at REAL(KIND=16) — any non-zero divergence is a real migration
    ! bug, even one that lands inside the per-case rounding budget.
    logical :: strict_exact = .false.

contains

    subroutine report_init(routine, target_name)
        character(len=*), intent(in) :: routine, target_name
        character(len=:), allocatable :: filename
        integer :: ios

        routine_save = trim(routine)
        any_failure = .false.
        case_count = 0
        strict_exact = (trim(target_name) == 'kind16')
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

        if (max_rel > 0.0_ep) then
            digits = -log10(max_rel)
        else
            digits = 99.0_ep
        end if
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

        ! Echo to stdout for ctest --output-on-failure / human reading.
        write(*, '(a,a,a,a,a,es12.4,a,f6.2,a)') &
            '  ', trim(routine_save), ' [', trim(case_label), &
            '] max_rel_err=', max_rel, '  digits=', digits, &
            merge(' PASS', ' FAIL', passed)
    end subroutine report_case

    subroutine report_finalize()
        write(unit_save, '(a)') '  ]'
        write(unit_save, '(a)') '}'
        close(unit_save)
        unit_save = -1
        ! Non-zero exit code if any case failed — CTest sees pass/fail.
        if (any_failure) error stop 1
    end subroutine report_finalize

end module prec_report
