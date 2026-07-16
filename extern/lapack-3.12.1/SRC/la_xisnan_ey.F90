!> \brief \b LA_XISNAN_EY provides LA_ISNAN overloads for extended
!> precision (KIND=10), complementing LA_XISNAN (single/double).
!
!  This is the kind10 half of the former LA_XISNAN_EP, split so each
!  extended target ships a single-precision module in its own prefixed
!  archive (E → libelapack) instead of a shared lapack_common. Because
!  the module now lives in the precision-specific archive — compiled only
!  when the kind10 target is built, which by definition needs KIND=10 —
!  no HAVE_REAL10 guard is required.
!
module LA_XISNAN_EY
   use LA_XISNAN

   interface LA_ISNAN
   module procedure EISNAN
   end interface

contains

   logical function EISNAN( x )
   use LA_CONSTANTS_EY, only: wp=>ep
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   EISNAN = ieee_is_nan(x)
#elif USE_ISNAN
   EISNAN = isnan(x)
#else
   EISNAN = ELAISNAN(x,x)

   contains
   logical function ELAISNAN( x, y )
   use LA_CONSTANTS_EY, only: wp=>ep
   real(wp) :: x, y
   ELAISNAN = ( x.ne.y )
   end function ELAISNAN
#endif
   end function EISNAN

end module LA_XISNAN_EY
