! Quad ↔ target-type converters for the kind8 baseline.
!
! kind8 is the unmigrated double-precision (REAL(KIND=8)) baseline:
! plain Netlib d/z routines compared against the quad reference.
! q2t_r casts down (loses precision relative to the quad reference);
! t2q_r promotes back up. The cast-down is the inherent floor of
! IEEE binary64 (~16 decimal digits of agreement) — a sanity-check
! the kind16 quad reference is faithful to that floor.

module target_conv
    implicit none
    private

    integer, parameter, public :: ep = 16
    integer, parameter, public :: tk = 8

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
