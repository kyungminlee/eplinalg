!> \brief \b LA_CONSTANTS_EY defines scaling constants for extended
!> precision (KIND=10, 80-bit x87), complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Naming conventions follow LA_CONSTANTS:
!    KIND=10:  E prefix (real), Y prefix (complex)
!
!  This is the kind10 half of the former LA_CONSTANTS_EP, split so each
!  extended target ships a single-precision module in its own prefixed
!  archive (E/Y → libelapack) instead of a shared lapack_common. Because
!  the module now lives in the precision-specific archive — compiled only
!  when the kind10 target is built, which by definition needs KIND=10 —
!  no HAVE_REAL10 guard is required.
!
module LA_CONSTANTS_EY
!  -- LAPACK auxiliary module (extended precision, KIND=10) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

! =====================================================================
!  Extended precision (KIND=10) — 80-bit x87 extended double
! =====================================================================
   integer, parameter :: ep = selected_real_kind(18, 4931)

!  Standard constants
   real(ep), parameter :: ezero = 0.0_ep
   real(ep), parameter :: ehalf = 0.5_ep
   real(ep), parameter :: eone = 1.0_ep
   real(ep), parameter :: etwo = 2.0_ep
   real(ep), parameter :: ethree = 3.0_ep
   real(ep), parameter :: efour = 4.0_ep
   real(ep), parameter :: eeight = 8.0_ep
   real(ep), parameter :: eten = 10.0_ep
   complex(ep), parameter :: yzero = ( 0.0_ep, 0.0_ep )
   complex(ep), parameter :: yhalf = ( 0.5_ep, 0.0_ep )
   complex(ep), parameter :: yone = ( 1.0_ep, 0.0_ep )
   character*1, parameter :: eprefix = 'E'
   character*1, parameter :: yprefix = 'Y'

!  Scaling constants
   real(ep), parameter :: eulp = epsilon(0._ep)
   real(ep), parameter :: eeps = eulp * 0.5_ep
   real(ep), parameter :: esafmin = real(radix(0._ep),ep)**max( &
      minexponent(0._ep)-1, &
      1-maxexponent(0._ep) &
   )
   real(ep), parameter :: esafmax = eone / esafmin
   real(ep), parameter :: esmlnum = esafmin / eulp
   real(ep), parameter :: ebignum = esafmax * eulp
   real(ep), parameter :: ertmin = sqrt(esmlnum)
   real(ep), parameter :: ertmax = sqrt(ebignum)

!  Blue's scaling constants
   real(ep), parameter :: etsml = real(radix(0._ep), ep)**ceiling( &
       (minexponent(0._ep) - 1) * 0.5_ep)
   real(ep), parameter :: etbig = real(radix(0._ep), ep)**floor( &
       (maxexponent(0._ep) - digits(0._ep) + 1) * 0.5_ep)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
   real(ep), parameter :: essml = real(radix(0._ep), ep)**( - floor( &
       (minexponent(0._ep) - digits(0._ep)) * 0.5_ep))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
   real(ep), parameter :: esbig = real(radix(0._ep), ep)**( - ceiling( &
       (maxexponent(0._ep) + digits(0._ep) - 1) * 0.5_ep))

end module LA_CONSTANTS_EY
