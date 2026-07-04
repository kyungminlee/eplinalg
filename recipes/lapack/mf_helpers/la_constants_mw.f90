!> \brief \b LA_CONSTANTS_MW defines scaling constants for the multifloats
!> double-double (real64x2) precision, complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Multifloats double-double (~106-bit mantissa) constants. Re-exports
!  named constants from the multifloats module under the migrator's
!  multifloats prefix (M for real, W for complex). The migrator
!  rewrites the USE-clause aliases in LAPACK source from
!  ``zero=>dzero`` (etc.) to ``zero=>mzero`` (etc.); the local LHS
!  alias is unchanged so the body of the routine compiles without
!  further substitution.
!
!  Naming conventions per migrator target:
!    KIND=8  (double):  D prefix (real), Z prefix (complex)
!    KIND=10 (extended): E prefix (real), Y prefix (complex)
!    KIND=16 (quad):    Q prefix (real), X prefix (complex)
!    multifloats:       M prefix (real), W prefix (complex)
!  Single-letter M/W are unused as routine prefixes anywhere in
!  upstream BLAS / LAPACK / ScaLAPACK, so the renamed namespace is
!  collision-free. (Earlier two-letter DD/ZZ prefixes collided with
!  ScaLAPACK's orphaned DDDOT wrapper, corrupting pddpotf2; earlier
!  T/V prefixes had local-name shadows in latrs / latbs / latps and
!  pXlahqr / pXlaqr1 that required per-file overrides.)
!
module LA_CONSTANTS_MW
   use multifloats, only: real64x2, cmplx64x2, &
                          DD_ZERO, DD_HALF, DD_ONE, DD_TWO, DD_EIGHT, &
                          DD_SAFMIN, DD_SAFMAX, &
                          DD_RTMIN, DD_RTMAX
   implicit none
   private
   public :: real64x2, cmplx64x2

   ! Local alias for the leading-limb double precision kind. Used to
   ! evaluate Blue's scaling constants from `radix/minexponent/maxexponent/digits`
   ! at compile time (see msbig/mssml/mtbig/mtsml below).
   integer, parameter :: dp = kind(1.0d0)

! =====================================================================
!  Multifloats double-double (~106-bit) constants
! =====================================================================

!  Standard real constants (M-prefixed to match the migrator rename of
!  the corresponding D-prefixed source symbols).
   type(real64x2), parameter, public :: mzero  = DD_ZERO
   type(real64x2), parameter, public :: mhalf  = DD_HALF
   type(real64x2), parameter, public :: mone   = DD_ONE
   type(real64x2), parameter, public :: mtwo   = DD_TWO
   ! Use named-component structure constructor (``limbs=...``) so that
   ! the compiler binds these initializers to the structure constructor
   ! of real64x2 rather than to the overloaded ``real64x2(...)``
   ! generic interface — the latter is a function call and is therefore
   ! illegal in a PARAMETER initializer.
   type(real64x2), parameter, public :: mthree = real64x2(limbs=[3.0d0, 0.0d0])
   type(real64x2), parameter, public :: mfour  = real64x2(limbs=[4.0d0, 0.0d0])
   type(real64x2), parameter, public :: meight = DD_EIGHT
   type(real64x2), parameter, public :: mten   = real64x2(limbs=[10.0d0, 0.0d0])

!  Complex constants (W-prefixed). Must use named-component structure
!  constructor syntax (``re=`` / ``im=``) so the compiler picks the
!  structure constructor and not the overloaded ``cmplx64x2`` interface
!  procedures (which are not allowed in PARAMETER initializers).
   type(cmplx64x2), parameter, public :: wzero = &
      cmplx64x2(re=DD_ZERO, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: whalf = &
      cmplx64x2(re=DD_HALF, im=DD_ZERO)
   type(cmplx64x2), parameter, public :: wone  = &
      cmplx64x2(re=DD_ONE,  im=DD_ZERO)

   character*1, parameter, public :: mprefix = 'M'
   character*1, parameter, public :: wprefix = 'W'

!  Scaling constants (mirror la_constants.f90 names with m prefix)
   type(real64x2), parameter, public :: msafmin = DD_SAFMIN
   type(real64x2), parameter, public :: msafmax = DD_SAFMAX
   type(real64x2), parameter, public :: mrtmin  = DD_RTMIN
   type(real64x2), parameter, public :: mrtmax  = DD_RTMAX

!  Blue's scaling constants — derived from the leading-limb double's
!  exponent range so `(ax * msbig)**2` stays representable for any
!  `ax` in [tiny(double), huge(double)]. NOT re-exported from
!  upstream multifloats's DD_TBIG/DD_SBIG/DD_TSML/DD_SSML — those are
!  defined as 1e±100 / 1e±50 (a SCALE-UP factor for sbig, opposite to
!  LAPACK's convention) and overflow `(ax * sbig)**2` once ax exceeds
!  ~1e+104, which dgesvj's protect-from-underflow up-scaling pass
!  routinely produces. The LAPACK formula:
!     tsml = base^ceiling((minexp-1)/2)
!     tbig = base^floor((maxexp-digits+1)/2)
!     ssml = base^(-floor((minexp-digits)/2))
!     sbig = base^(-ceiling((maxexp+digits-1)/2))
!  applied at the leading-limb kind (`dp` = double) gives full-range
!  safety on a DD-on-double type (the lo limb is below the hi limb's
!  ULP, so the hi limb's exponent dictates overflow boundaries).
   type(real64x2), parameter, public :: mtsml = real64x2(limbs=[ &
      real(radix(0.0_dp), dp)**ceiling( &
         (minexponent(0.0_dp) - 1) * 0.5_dp), 0.0_dp])
   type(real64x2), parameter, public :: mtbig = real64x2(limbs=[ &
      real(radix(0.0_dp), dp)**floor( &
         (maxexponent(0.0_dp) - digits(0.0_dp) + 1) * 0.5_dp), 0.0_dp])
   type(real64x2), parameter, public :: mssml = real64x2(limbs=[ &
      real(radix(0.0_dp), dp)**(-floor( &
         (minexponent(0.0_dp) - digits(0.0_dp)) * 0.5_dp)), 0.0_dp])
   type(real64x2), parameter, public :: msbig = real64x2(limbs=[ &
      real(radix(0.0_dp), dp)**(-ceiling( &
         (maxexponent(0.0_dp) + digits(0.0_dp) - 1) * 0.5_dp)), 0.0_dp])

end module LA_CONSTANTS_MW
