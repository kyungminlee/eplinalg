! Quad ↔ target-type converters for the kind4 baseline.
!
! kind4 is the unmigrated single-precision (REAL(KIND=4)) baseline:
! plain Netlib s/c routines compared against the quad reference.
! q2t_r casts down to single (~7 decimal digits); t2q_r promotes back
! up. The kind4 column in tests/RESULT.md is both a precision floor
! and an independent sanity-check on the quad reference — it touches
! none of the -freal-8-real-16 promotion machinery, so a kind4 row
! that disagrees by far more than ~7 digits flags an issue with the
! reference rather than with the migrated target.

module target_conv
    implicit none
    private

    integer, parameter, public :: ep = 16
    integer, parameter, public :: tk = 4

    public :: q2t_r, t2q_r, q2t_c, t2q_c

contains

    elemental function q2t_r(x) result(y)
        real(ep), intent(in) :: x
        real(tk) :: y
        y = real(x, tk)
    end function

    elemental function t2q_r(x) result(y)
        real(tk), intent(in) :: x
        real(ep) :: y
        y = real(x, ep)
    end function

    elemental function q2t_c(z) result(w)
        complex(ep), intent(in) :: z
        complex(tk) :: w
        w = cmplx(z, kind=tk)
    end function

    elemental function t2q_c(z) result(w)
        complex(tk), intent(in) :: z
        complex(ep) :: w
        w = cmplx(z, kind=ep)
    end function

end module target_conv
