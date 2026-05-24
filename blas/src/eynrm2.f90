!> \brief \b EYNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL(KIND=10) FUNCTION EYNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX(KIND=10) X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> EYNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    EYNRM2 := sqrt( x**H*x )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX(KIND=10) array, dimension (N)
!>         complex vector with N elements
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup nrm2
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
!
!  Override of upstream LAPACK 3.12.1's BLAS/SRC/eynrm2.f90. See the
!  matching enrm2.f90 override for the full rationale: Blue's scaling
!  constants are renamed `tsml/tbig/ssml/sbig` -> `btsml/btbig/bssml/bsbig`
!  (avoids the multifloats-target nuke-rename that imports broken
!  `DD_T*BIG / DD_S*ML` values), and they are SAVE locals init'd at
!  first call rather than PARAMETERs (multifloats `**` operator is
!  user-defined and not allowed in initializer expressions).
function EYNRM2( n, x, incx )
   integer, parameter :: wp = 10
   real(wp) :: EYNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
!  ..
!  .. Blue's scaling constants + maxN (SAVE; init at first call below).
!     maxN is a SAVE rather than PARAMETER because the multifloats target
!     overloads `huge` for `real64x2`, and gfortran rejects user-defined
!     functions in PARAMETER initializer expressions.
   real(wp), save :: maxN
   real(wp), save :: btsml, btbig, bssml, bsbig
   logical, save :: blue_initialized = .false.
   real(wp) :: r_radix, half_wp
   integer  :: minexp_w, maxexp_w, digits_w
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  First-call init of Blue's scaling constants (see enrm2 override).
!
   if (.not. blue_initialized) then
      r_radix  = real(radix(0._wp), wp)
      half_wp  = 0.5_wp
      minexp_w = minexponent(0._wp)
      maxexp_w = maxexponent(0._wp)
      digits_w = digits(0._wp)
      btsml = r_radix**ceiling((minexp_w - 1) * half_wp)
      btbig = r_radix**floor((maxexp_w - digits_w + 1) * half_wp)
      bssml = r_radix**(-floor((minexp_w - digits_w) * half_wp))
      bsbig = r_radix**(-ceiling((maxexp_w + digits_w - 1) * half_wp))
      maxN = huge(0._wp)
      blue_initialized = .true.
   end if
!
!  Quick return if possible
!
   EYNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     btbig -- values bigger than this are scaled down by bsbig
!     btsml -- values smaller than this are scaled up by bssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(REAL(x(ix), KIND=10))
      if (ax > btbig) then
         abig = abig + (ax*bsbig)**2
         notbig = .false.
      else if (ax < btsml) then
         if (notbig) asml = asml + (ax*bssml)**2
      else
         amed = amed + ax**2
      end if
      ax = abs(aimag(x(ix)))
      if (ax > btbig) then
         abig = abig + (ax*bsbig)**2
         notbig = .false.
      else if (ax < btsml) then
         if (notbig) asml = asml + (ax*bssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*bsbig)*bsbig
      end if
      scl = one / bsbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / bssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / bssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   EYNRM2 = scl*sqrt( sumsq )
   return
end function
