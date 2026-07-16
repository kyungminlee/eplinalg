! Quad ↔ target-type converters for the kind16 build.
!
! kind16's working type is REAL(KIND=16) — the same kind as the
! reference precision `ep`, so the converters are identities. The
! module exists so the templated test wrappers can call q2t_r / t2q_r /
! q2t_c / t2q_c uniformly across all targets without per-target
! conditionals in the body.

module target_conv
    implicit none
    private

    integer, parameter, public :: ep = 16

    public :: q2t_r, t2q_r, q2t_c, t2q_c

contains

    elemental function q2t_r(x) result(y)
        real(ep), intent(in) :: x
        real(ep) :: y
        y = x
    end function

    elemental function t2q_r(x) result(y)
        real(ep), intent(in) :: x
        real(ep) :: y
        y = x
    end function

    elemental function q2t_c(z) result(w)
        complex(ep), intent(in) :: z
        complex(ep) :: w
        w = z
    end function

    elemental function t2q_c(z) result(w)
        complex(ep), intent(in) :: z
        complex(ep) :: w
        w = z
    end function

end module target_conv
