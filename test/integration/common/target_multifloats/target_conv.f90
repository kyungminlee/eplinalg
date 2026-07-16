! Quad ↔ target-type converters for the multifloats build.
!
! Multifloats uses a double-double (TYPE(real64x2)) representation of
! ~32 decimal digits. q2t_r splits each quad value into the high
! (double approximation) and low (double remainder) limbs:
!     hi = round-to-double(quad)
!     lo = round-to-double(quad - hi)
! `lo` fits in a double because |lo| <= ulp(hi)/2 < 2^-52 * |hi|.
! t2q_r recombines the two limbs back into a single quad value.

module target_conv
    use multifloats, only: real64x2, cmplx64x2
    implicit none
    private

    integer, parameter, public :: ep = 16
    integer, parameter, public :: dp = 8

    public :: q2t_r, t2q_r, q2t_c, t2q_c

contains

    elemental function q2t_r(x) result(y)
        real(ep), intent(in) :: x
        type(real64x2) :: y
        real(dp) :: hi
        hi = real(x, dp)
        y%limbs(1) = hi
        y%limbs(2) = real(x - real(hi, ep), dp)
    end function

    elemental function t2q_r(x) result(y)
        type(real64x2), intent(in) :: x
        real(ep) :: y
        y = real(x%limbs(1), ep) + real(x%limbs(2), ep)
    end function

    elemental function q2t_c(z) result(w)
        complex(ep), intent(in) :: z
        type(cmplx64x2) :: w
        w%re = q2t_r(real(z, ep))
        w%im = q2t_r(aimag(z))
    end function

    elemental function t2q_c(z) result(w)
        type(cmplx64x2), intent(in) :: z
        complex(ep) :: w
        w = cmplx(t2q_r(z%re), t2q_r(z%im), ep)
    end function

end module target_conv
