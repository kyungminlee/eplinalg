!> \brief \b LA_CONSTANTS_QX defines scaling constants for quadruple
!> precision (KIND=16, 128-bit IEEE binary128), complementing LA_CONSTANTS.
!
!  =========== DOCUMENTATION ===========
!
!  Naming conventions follow LA_CONSTANTS:
!    KIND=16:  Q prefix (real), X prefix (complex)
!
!  This is the kind16 half of the former LA_CONSTANTS_EP, split so each
!  extended target ships a single-precision module in its own prefixed
!  archive (Q/X → libqlapack) instead of a shared lapack_common. Because
!  the module now lives in the precision-specific archive — compiled only
!  when the kind16 target is built, which by definition needs KIND=16 —
!  no HAVE_REAL16 guard is required.
!
module LA_CONSTANTS_QX
!  -- LAPACK auxiliary module (quadruple precision, KIND=16) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

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

end module LA_CONSTANTS_QX
