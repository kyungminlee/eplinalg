!> \brief \b LA_XISNAN_QX provides LA_ISNAN overloads for quadruple
!> precision (KIND=16), complementing LA_XISNAN (single/double).
!
!  This is the kind16 half of the former LA_XISNAN_EP, split so each
!  extended target ships a single-precision module in its own prefixed
!  archive (Q → libqlapack) instead of a shared lapack_common. Because
!  the module now lives in the precision-specific archive — compiled only
!  when the kind16 target is built, which by definition needs KIND=16 —
!  no HAVE_REAL16 guard is required.
!
module LA_XISNAN_QX
   use LA_XISNAN

   interface LA_ISNAN
   module procedure QISNAN
   end interface

contains

   logical function QISNAN( x )
   use LA_CONSTANTS_QX, only: wp=>qp
#ifdef USE_IEEE_INTRINSIC
   use, intrinsic :: ieee_arithmetic
#elif USE_ISNAN
   intrinsic :: isnan
#endif
   real(wp) :: x
#ifdef USE_IEEE_INTRINSIC
   QISNAN = ieee_is_nan(x)
#elif USE_ISNAN
   QISNAN = isnan(x)
#else
   QISNAN = QLAISNAN(x,x)

   contains
   logical function QLAISNAN( x, y )
   use LA_CONSTANTS_QX, only: wp=>qp
   real(wp) :: x, y
   QLAISNAN = ( x.ne.y )
   end function QLAISNAN
#endif
   end function QISNAN

end module LA_XISNAN_QX
