      SUBROUTINE YLANV2( A, B, C, D, RT1, RT2, CS, SN )
      IMPLICIT NONE
*
*  -- ScaLAPACK routine (version 1.7) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     May 28, 1999
*
*     .. Scalar Arguments ..
      REAL(KIND=10)   CS
      COMPLEX(KIND=10)         A, B, C, D, RT1, RT2, SN
*     ..
*
*  Purpose
*  =======
*
*  YLANV2 computes the Schur factorization of a complex 2-by-2
*  nonhermitian matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
*
*  Arguments
*  =========
*
*  A       (input/output) COMPLEX(KIND=10)
*  B       (input/output) COMPLEX(KIND=10)
*  C       (input/output) COMPLEX(KIND=10)
*  D       (input/output) COMPLEX(KIND=10)
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1     (output) COMPLEX(KIND=10)
*  RT2     (output) COMPLEX(KIND=10)
*          The two eigenvalues.
*
*  CS      (output) REAL(KIND=10)
*  SN      (output) COMPLEX(KIND=10)
*          Parameters of the rotation matrix.
*
*  Further Details
*  ===============
*
*  Implemented by Mark R. Fahey, May 28, 1999
*
*  =====================================================================
*
*     .. Parameters ..
      REAL(KIND=10)   RZERO, HALF, RONE
      PARAMETER          ( RZERO = 0.0E0_10, HALF = 0.5E0_10,           
     +         RONE = 1.0E0_10 )
      COMPLEX(KIND=10)         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0E0_10, 0.0E0_10 ),               
     +     ONE = ( 1.0E0_10, 0.0E0_10 ) )
*     ..
*     .. Local Scalars ..
      COMPLEX(KIND=10)         AA, BB, DD, T, TEMP, TEMP2, U, X, Y
      REAL(KIND=10)   ZR, ZI
*     ..
*     .. External Subroutines ..
      EXTERNAL           YLARTG, ELADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Initialize CS and SN
*
      CS = RONE
      SN = ZERO
*
      IF( C.EQ.ZERO ) THEN
         GO TO 10
*
      ELSE IF( B.EQ.ZERO ) THEN
*
*        Swap rows and columns
*
         CS = RZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( ( A-D ).EQ.ZERO ) THEN
         TEMP = SQRT( B*C )
         A = A + TEMP
         D = D - TEMP
         IF( ( B+C ).EQ.ZERO ) THEN
            CS = SQRT( HALF )
            SN = CMPLX( RZERO, RONE , KIND=10)*CS
         ELSE
            TEMP = SQRT( B+C )
            CALL ELADIV( REAL( SQRT( B ) , KIND=10), AIMAG( SQRT( B ) ),
     +                    REAL( TEMP , KIND=10), AIMAG( TEMP ), ZR, ZI )
            TEMP2 = CMPLX( ZR, ZI , KIND=10)
            CS = REAL( TEMP2 , KIND=10)
            CALL ELADIV( REAL( SQRT( C ) , KIND=10), AIMAG( SQRT( C ) ),
     +                    REAL( TEMP , KIND=10), AIMAG( TEMP ), ZR, ZI )
            SN = CMPLX( ZR, ZI , KIND=10)
         END IF
         B = B - C
         C = ZERO
         GO TO 10
      ELSE
*
*        Compute eigenvalue closest to D
*
         T = D
         U = B*C
         X = HALF*( A-T )
         Y = SQRT( X*X+U )
         IF( REAL( X , KIND=10)*REAL( Y , KIND=10)+AIMAG( X )*AIMAG( Y 
     +).LT.RZERO )       Y = -Y
         CALL ELADIV( REAL( U , KIND=10), AIMAG( U ),                 
     +REAL( X+Y , KIND=10), AIMAG( X+Y ), ZR, ZI )
         T = T - CMPLX( ZR, ZI , KIND=10)
*
*        Do one QR step with exact shift T - resulting 2 x 2 in
*        triangular form.
*
         CALL YLARTG( A-T, C, CS, SN, AA )
*
         D = D - T
         BB = CS*B + SN*D
         DD = -CONJG( SN )*B + CS*D
*
         A = AA*CS + BB*CONJG( SN ) + T
         B = -AA*SN + BB*CS
         C = ZERO
         D = T
*
      END IF
*
   10 CONTINUE
*
*     Store eigenvalues in RT1 and RT2.
*
      RT1 = A
      RT2 = D
      RETURN
*
*     End of YLANV2
*
      END
