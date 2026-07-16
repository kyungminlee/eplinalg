! Quad ↔ target-type converters for the kind10 build.
!
! kind10's working type is REAL(KIND=10) — narrower than the reference
! quad. q2t_r casts down (loses precision); t2q_r promotes back up
! (preserves the now-rounded value). The cast-down is the inherent
! precision floor of this target; comparison against the quad reference
! measures how closely kind10 reproduces the quad value.

module target_conv
    implicit none
    private

    integer, parameter, public :: ep = 16
    integer, parameter, public :: tk = 10

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
