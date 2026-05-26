!> \brief \b LA_XISNAN_EP provides LA_ISNAN overloads for extended and
!> quadruple precision, complementing LA_XISNAN (single/double).
!
!  Compile with -DHAVE_REAL10 to include the KIND=10 overload.
!
module LA_XISNAN_EP
   use LA_XISNAN
   interface LA_ISNAN

#ifdef HAVE_REAL10
   module procedure EISNAN
#endif
#ifdef HAVE_REAL16
   module procedure QISNAN
#endif

   end interface

contains

#ifdef HAVE_REAL10
   logical function EISNAN( x )
   use LA_CONSTANTS_EP, only: wp=>ep
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
   use LA_CONSTANTS_EP, only: wp=>ep
   real(wp) :: x, y
   ELAISNAN = ( x.ne.y )
   end function ELAISNAN
#endif
   end function EISNAN
#endif

#ifdef HAVE_REAL16
   logical function QISNAN( x )
   use LA_CONSTANTS_EP, only: wp=>qp
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
   use LA_CONSTANTS_EP, only: wp=>qp
   real(wp) :: x, y
   QLAISNAN = ( x.ne.y )
   end function QLAISNAN
#endif
   end function QISNAN
#endif

end module LA_XISNAN_EP
