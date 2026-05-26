*> \brief \b ELARF1L applies an elementary reflector to a general rectangular
*              matrix assuming v(lastv) = 1 where lastv is the last non-zero
*              element
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ELARF1L + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/elarf1l.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/elarf1l.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/elarf1l.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ELARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          SIDE
*       INTEGER            INCV, LDC, M, N
*       REAL(KIND=10)   TAU
*       ..
*       .. Array Arguments ..
*       REAL(KIND=10)   C( LDC, * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ELARF1L applies a real elementary reflector H to a real m by n matrix
*> C, from either the left or the right. H is represented in the form
*>
*>       H = I - tau * v * v**T
*>
*> where tau is a real scalar and v is a real vector.
*>
*> If tau = 0, then H is taken to be the unit matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': form  H * C
*>          = 'R': form  C * H
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is REAL(KIND=10) array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is REAL(KIND=10)
*>          The value tau in the representation of H.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is REAL(KIND=10) array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*>          or C * H if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL(KIND=10) array, dimension
*>                         (N) if SIDE = 'L'
*>                      or (M) if SIDE = 'R'
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup larf
*
*  =====================================================================
      SUBROUTINE ELARF1L( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      REAL(KIND=10)   TAU
*     ..
*     .. Array Arguments ..
      REAL(KIND=10)   C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   ONE, ZERO
      PARAMETER          ( ONE = 1.0E0_10, ZERO = 0.0E0_10 )
*     ..
*     .. Local Scalars ..
      LOGICAL            APPLYLEFT
      INTEGER            I, FIRSTV, LASTV, LASTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           EAXPY, EGEMV, EGER, ESCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAELR, ILAELC
      EXTERNAL           LSAME, ILAELR, ILAELC
*     ..
*     .. Executable Statements ..
*
      APPLYLEFT = LSAME( SIDE, 'L' )
      FIRSTV = 1
      LASTC = 0
      IF( TAU.NE.ZERO ) THEN
!     Set up variables for scanning V.  LASTV begins pointing to the end
!     of V.
         IF( APPLYLEFT ) THEN
            LASTV = M
         ELSE
            LASTV = N
         END IF
         I = 1
!     Look for the last non-zero row in V.
         DO WHILE( LASTV.GT.FIRSTV .AND. V( I ).EQ.ZERO )
            FIRSTV = FIRSTV + 1
            I = I + INCV
         END DO
         IF( APPLYLEFT ) THEN
!     Scan for the last non-zero column in C(1:lastv,:).
            LASTC = ILAELC(LASTV, N, C, LDC)
         ELSE
!     Scan for the last non-zero row in C(:,1:lastv).
            LASTC = ILAELR(M, LASTV, C, LDC)
         END IF
      END IF
      IF( LASTC.EQ.0 ) THEN
         RETURN
      END IF
      IF( APPLYLEFT ) THEN
*
*        Form  H * C
*
         IF( LASTV.GT.0 ) THEN
            ! Check if m = 1.E0_10 This means v = 1, So we just need to compute
            ! C := HC = (1-\tau)C.
            IF( LASTV.EQ.FIRSTV ) THEN
               CALL ESCAL(LASTC, ONE - TAU, C( FIRSTV, 1), LDC)
            ELSE
*
*              w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
*
               ! w(1:lastc,1) := C(1:lastv-1,1:lastc)**T * v(1:lastv-1,1)
               CALL EGEMV( 'Transpose', LASTV-FIRSTV, LASTC, ONE,       
     +               C(FIRSTV,1), LDC, V(I), INCV, ZERO,                
     +       WORK, 1)
               ! w(1:lastc,1) += C(lastv,1:lastc)**T * v(lastv,1) = C(lastv,1:lastc)**T
               CALL EAXPY(LASTC, ONE, C(LASTV,1), LDC, WORK, 1)
*
*              C(1:lastv,1:lastc) := C(...) - tau * v(1:lastv,1) * w(1:lastc,1)**T
*
               ! C(lastv, 1:lastc)   := C(...) - tau * v(lastv,1) * w(1:lastc,1)**T
               !                      = C(...) - tau * w(1:lastc,1)**T
               CALL EAXPY(LASTC, -TAU, WORK, 1, C(LASTV,1), LDC)
               ! C(1:lastv-1,1:lastc) := C(...) - tau * v(1:lastv-1,1)*w(1:lastc,1)**T
               CALL EGER(LASTV-FIRSTV, LASTC, -TAU, V(I), INCV,         
     +           WORK, 1, C(FIRSTV,1), LDC)
            END IF
         END IF
      ELSE
*
*        Form  C * H
*
         IF( LASTV.GT.0 ) THEN
            ! Check if n = 1.E0_10 This means v = 1, so we just need to compute
            ! C := CH = C(1-\tau).
            IF( LASTV.EQ.FIRSTV ) THEN
               CALL ESCAL(LASTC, ONE - TAU, C, 1)
            ELSE
*
*              w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
*
               ! w(1:lastc,1) := C(1:lastc,1:lastv-1) * v(1:lastv-1,1)
               CALL EGEMV( 'No transpose', LASTC, LASTV-FIRSTV,         
     +    ONE, C(1,FIRSTV), LDC, V(I), INCV, ZERO, WORK, 1 )
               ! w(1:lastc,1) += C(1:lastc,lastv) * v(lastv,1) = C(1:lastc,lastv)
               CALL EAXPY(LASTC, ONE, C(1,LASTV), 1, WORK, 1)
*
*              C(1:lastc,1:lastv) := C(...) - tau * w(1:lastc,1) * v(1:lastv,1)**T
*
               ! C(1:lastc,lastv)     := C(...) - tau * w(1:lastc,1) * v(lastv,1)**T
               !                       = C(...) - tau * w(1:lastc,1)
               CALL EAXPY(LASTC, -TAU, WORK, 1, C(1,LASTV), 1)
               ! C(1:lastc,1:lastv-1) := C(...) - tau * w(1:lastc,1) * v(1:lastv-1)**T
               CALL EGER( LASTC, LASTV-FIRSTV, -TAU, WORK, 1, V(I),     
     +                  INCV, C(1,FIRSTV), LDC )
            END IF
         END IF
      END IF
      RETURN
*
*     End of ELARF1L
*
      END
