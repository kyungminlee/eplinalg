!> \brief \b LA_CONSTANTS_EP defines scaling constants for extended and
!> quadruple precision, complementing LA_CONSTANTS (single/double).
!
!  =========== DOCUMENTATION ===========
!
!  Extended-precision (KIND=10, 80-bit) constants are only available
!  when the compiler and target architecture support them.  Compile
!  with -DHAVE_REAL10 to enable.  GFortran on x86/x86_64 supports
!  this; Intel Fortran (ifort/ifx) on x86 does not.
!
!  Quadruple-precision (KIND=16, 128-bit) constants are always
!  available on compilers that support REAL(16) (GFortran, Intel
!  Fortran, NAG, etc.).
!
!  Naming conventions follow LA_CONSTANTS:
!    KIND=10:  E prefix (real), Y prefix (complex)
!    KIND=16:  Q prefix (real), X prefix (complex)
!
module LA_CONSTANTS_EP
!  -- LAPACK auxiliary module (extended / quad precision) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

#ifdef HAVE_REAL10
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

#endif

#ifdef HAVE_REAL16
! =====================================================================
!  Quadruple precision (KIND=16) — 128-bit IEEE binary128
! =====================================================================
   integer, parameter :: qp = selected_real_kind(33)

!  Standard constants
   real(qp), parameter :: qzero = 0.0_qp
   real(qp), parameter :: qhalf = 0.5_qp
   real(qp), parameter :: qone = 1.0_qp
   real(qp), parameter :: qtwo = 2.0_qp
   real(qp), parameter :: qthree = 3.0_qp
   real(qp), parameter :: qfour = 4.0_qp
   real(qp), parameter :: qeight = 8.0_qp
   real(qp), parameter :: qten = 10.0_qp
   complex(qp), parameter :: xzero = ( 0.0_qp, 0.0_qp )
   complex(qp), parameter :: xhalf = ( 0.5_qp, 0.0_qp )
   complex(qp), parameter :: xone = ( 1.0_qp, 0.0_qp )
   character*1, parameter :: qprefix = 'Q'
   character*1, parameter :: xprefix = 'X'

!  Scaling constants
   real(qp), parameter :: qulp = epsilon(0._qp)
   real(qp), parameter :: qeps = qulp * 0.5_qp
   real(qp), parameter :: qsafmin = real(radix(0._qp),qp)**max( &
      minexponent(0._qp)-1, &
      1-maxexponent(0._qp) &
   )
   real(qp), parameter :: qsafmax = qone / qsafmin
   real(qp), parameter :: qsmlnum = qsafmin / qulp
   real(qp), parameter :: qbignum = qsafmax * qulp
   real(qp), parameter :: qrtmin = sqrt(qsmlnum)
   real(qp), parameter :: qrtmax = sqrt(qbignum)

!  Blue's scaling constants
   real(qp), parameter :: qtsml = real(radix(0._qp), qp)**ceiling( &
       (minexponent(0._qp) - 1) * 0.5_qp)
   real(qp), parameter :: qtbig = real(radix(0._qp), qp)**floor( &
       (maxexponent(0._qp) - digits(0._qp) + 1) * 0.5_qp)
!  ssml >= 1/s, where s was defined in https://doi.org/10.1145/355769.355771
!  The correction was added in https://doi.org/10.1145/3061665 to scale denormalized numbers correctly
   real(qp), parameter :: qssml = real(radix(0._qp), qp)**( - floor( &
       (minexponent(0._qp) - digits(0._qp)) * 0.5_qp))
!  sbig = 1/S, where S was defined in https://doi.org/10.1145/355769.355771
   real(qp), parameter :: qsbig = real(radix(0._qp), qp)**( - ceiling( &
       (maxexponent(0._qp) + digits(0._qp) - 1) * 0.5_qp))

#endif

end module LA_CONSTANTS_EP
