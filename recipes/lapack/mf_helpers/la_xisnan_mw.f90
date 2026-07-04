!> \brief \b LA_XISNAN_MW — multifloats analogue of LA_XISNAN.
!
!  =========== DOCUMENTATION ===========
!
!  Provides ``LA_ISNAN`` for ``type(real64x2)`` arguments. A
!  double-double value is NaN iff its first (high-order) limb is NaN —
!  the second limb is normally zero or carries low-order bits, and any
!  arithmetic that produces a NaN propagates it through the high limb.
!
module LA_XISNAN_MW
   use multifloats, only: real64x2
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private
   public :: LA_ISNAN

   interface LA_ISNAN
      module procedure WISNAN
   end interface

contains

   logical function WISNAN(x)
      type(real64x2), intent(in) :: x
      WISNAN = ieee_is_nan(x%limbs(1))
   end function WISNAN

end module LA_XISNAN_MW
