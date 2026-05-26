C
C  This file is part of MUMPS 5.8.2, released
C  on Mon Jan 12 15:17:08 UTC 2026
C
C
C  Copyright 1991-2026 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      SUBROUTINE EMUMPS_UPDATEDETER(PIV, DETER, NEXP)
      IMPLICIT NONE
      REAL(KIND=10), intent(in) :: PIV
      REAL(KIND=10), intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE EMUMPS_UPDATEDETER
      SUBROUTINE EMUMPS_UPDATEDETER_SCALING(PIV, DETER, NEXP)
      IMPLICIT NONE
      REAL(KIND=10), intent(in) :: PIV
      REAL(KIND=10), intent(inout) :: DETER
      INTEGER, intent(inout) :: NEXP
      DETER=DETER*fraction(PIV)
      NEXP=NEXP+exponent(PIV)+exponent(DETER)
      DETER=fraction(DETER)
      RETURN
      END SUBROUTINE EMUMPS_UPDATEDETER_SCALING
      SUBROUTINE EMUMPS_GETDETER2D(BLOCK_SIZE,IPIV,                     
     +  MYROW, MYCOL, NPROW, NPCOL,                       A, LOCAL_M, 
     +LOCAL_N, N, MYID,                       DETER,NEXP,SYM)
      IMPLICIT NONE
      INTEGER, intent (in)    :: SYM
      INTEGER, intent (inout) :: NEXP
      REAL(KIND=10), intent (inout) :: DETER
      INTEGER, intent (in)    :: BLOCK_SIZE, NPROW, NPCOL,
     &                           LOCAL_M, LOCAL_N, N
      INTEGER, intent (in)    :: MYROW, MYCOL, MYID, IPIV(LOCAL_M)
      REAL(KIND=10), intent(in) :: A(*)
      INTEGER  I,IMX,DI,NBLOCK,IBLOCK,ILOC,JLOC,
     &         ROW_PROC,COL_PROC, K
      DI = LOCAL_M + 1 
      NBLOCK = ( N - 1 ) / BLOCK_SIZE
      DO IBLOCK = 0, NBLOCK
        ROW_PROC = mod( IBLOCK, NPROW ) 
        IF ( MYROW.EQ.ROW_PROC ) THEN
          COL_PROC = mod( IBLOCK, NPCOL )
          IF ( MYCOL.EQ.COL_PROC ) THEN
            ILOC = ( IBLOCK / NPROW ) * BLOCK_SIZE
            JLOC = ( IBLOCK / NPCOL ) * BLOCK_SIZE
            I   =   ILOC + JLOC *  LOCAL_M + 1
            IMX = min(ILOC+BLOCK_SIZE,LOCAL_M)
     &            + (min(JLOC+BLOCK_SIZE,LOCAL_N)-1)*LOCAL_M
     &            + 1
            K=1
            DO WHILE ( I .LT. IMX )
              CALL EMUMPS_UPDATEDETER(A(I),DETER,NEXP)
              IF (SYM.EQ.1) THEN
               CALL EMUMPS_UPDATEDETER(A(I),DETER,NEXP)
              ENDIF
              IF (SYM.NE.1) THEN 
                IF (IPIV(ILOC+K) .NE. IBLOCK*BLOCK_SIZE+K) THEN
                  DETER = -DETER
                ENDIF
              ENDIF
              K = K + 1
              I = I + DI
            END DO
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EMUMPS_GETDETER2D
      SUBROUTINE EMUMPS_DETER_REDUCTION(            COMM, DETER_IN, 
     +NEXP_IN,            DETER_OUT, NEXP_OUT, NPROCS)
      IMPLICIT NONE
      INTEGER, intent(in) :: COMM, NPROCS
      REAL(KIND=10), intent(in) :: DETER_IN
      INTEGER,intent(in) :: NEXP_IN
      REAL(KIND=10),intent(out):: DETER_OUT
      INTEGER,intent(out):: NEXP_OUT
      INTEGER            :: IERR_MPI
      EXTERNAL EMUMPS_DETERREDUCE_FUNC
      INTEGER TWO_SCALARS_TYPE, DETERREDUCE_OP
      REAL(KIND=10) :: INV(2)
      REAL(KIND=10) :: OUTV(2)
      INCLUDE 'mpif.h'
      IF (NPROCS .EQ. 1) THEN
        DETER_OUT = DETER_IN
        NEXP_OUT  = NEXP_IN
        RETURN
      ENDIF
      CALL MPI_TYPE_CONTIGUOUS(2, MPI_LONG_DOUBLE,
     &                         TWO_SCALARS_TYPE,
     &                         IERR_MPI)
      CALL MPI_TYPE_COMMIT(TWO_SCALARS_TYPE, IERR_MPI)
      CALL MPI_OP_CREATE(EMUMPS_DETERREDUCE_FUNC,                      
     +.TRUE.,                                 DETERREDUCE_OP,           
     +              IERR_MPI)
      INV(1)=DETER_IN
      INV(2)=REAL(NEXP_IN, KIND=10)
      CALL MPI_ALLREDUCE( INV, OUTV, 1, TWO_SCALARS_TYPE,
     &                    DETERREDUCE_OP, COMM, IERR_MPI)
      CALL MPI_OP_FREE(DETERREDUCE_OP, IERR_MPI)
      CALL MPI_TYPE_FREE(TWO_SCALARS_TYPE, IERR_MPI)
      DETER_OUT = OUTV(1)
      NEXP_OUT  = int(OUTV(2))
      RETURN
      END SUBROUTINE EMUMPS_DETER_REDUCTION
      SUBROUTINE EMUMPS_DETERREDUCE_FUNC(INV, INOUTV, NEL, DATATYPE)
      IMPLICIT NONE
#if defined(WORKAROUNDINTELILP64MPI2INTEGER) || defined(WORKAROUNDILP64MPICUSTOMREDUCE)
      INTEGER(4), INTENT(IN)    :: NEL, DATATYPE
#else
      INTEGER, INTENT(IN)    :: NEL, DATATYPE
#endif
      REAL(KIND=10), INTENT(IN)    :: INV    ( 2 * NEL )
      REAL(KIND=10), INTENT(INOUT) :: INOUTV ( 2 * NEL )
      INTEGER I, TMPEXPIN, TMPEXPINOUT
      DO I = 1, NEL
        TMPEXPIN    = int(INV   (I*2))
        TMPEXPINOUT = int(INOUTV(I*2))
        CALL EMUMPS_UPDATEDETER(INV(I*2-1),                            
     +INOUTV(I*2-1),                            TMPEXPINOUT)   
        TMPEXPINOUT = TMPEXPINOUT + TMPEXPIN
        INOUTV(I*2) = REAL(TMPEXPINOUT, KIND=10)
      ENDDO
      RETURN
      END SUBROUTINE EMUMPS_DETERREDUCE_FUNC
      SUBROUTINE EMUMPS_DETER_SQUARE(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      REAL(KIND=10), intent (inout) :: DETER
      DETER=DETER*DETER
      NEXP=NEXP+NEXP
      RETURN
      END SUBROUTINE EMUMPS_DETER_SQUARE
      SUBROUTINE EMUMPS_DETER_SCALING_INVERSE(DETER, NEXP)
      IMPLICIT NONE
      INTEGER, intent (inout) :: NEXP
      REAL(KIND=10), intent (inout) :: DETER
      DETER=1.0E0_10/DETER
      NEXP=-NEXP
      RETURN
      END SUBROUTINE EMUMPS_DETER_SCALING_INVERSE
      SUBROUTINE EMUMPS_DETER_SIGN_PERM(DETER, N, PERM)
      IMPLICIT NONE
      REAL(KIND=10), intent(inout) :: DETER
      INTEGER, intent(in)    :: N
      INTEGER, intent(inout) :: PERM(N)
      INTEGER I, J, K
      K = 0
      DO I = 1, N
        IF (PERM(I) .LT. 0) THEN
          PERM(I)=-PERM(I)
        ELSE 
          J = PERM(I)
          DO WHILE (J.NE.I)
            PERM(J)=-PERM(J)
            K = K + 1
            J = -PERM(J) 
          ENDDO
        ENDIF
      ENDDO
      IF (mod(K,2).EQ.1) THEN
        DETER = -DETER
      ENDIF
      RETURN
      END SUBROUTINE EMUMPS_DETER_SIGN_PERM
      SUBROUTINE EMUMPS_PAR_ROOT_MINMAX_PIV_UPD (                       
     +BLOCK_SIZE,IPIV,                       MYROW, MYCOL, NPROW, NPCOL,
     +                       A, LOCAL_M, LOCAL_N, N, MYID,              
     +         DKEEP, KEEP, SYM)
       USE EMUMPS_FAC_FRONT_AUX_M,        ONLY : 
     +EMUMPS_UPDATE_MINMAX_PIVOT 
      IMPLICIT NONE
      INTEGER, intent (in)    :: BLOCK_SIZE, NPROW, NPCOL,
     &                           LOCAL_M, LOCAL_N, N, SYM
      INTEGER, intent (in)    :: MYROW, MYCOL, MYID, IPIV(LOCAL_M)
      REAL(KIND=10), intent(in) :: A(*)
      REAL(KIND=10), INTENT(INOUT) :: DKEEP(230)
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER  I,IMX,DI,NBLOCK,IBLOCK,ILOC,JLOC,
     &         ROW_PROC,COL_PROC, K
      REAL(KIND=10) :: ABSPIVOT
      DI = LOCAL_M + 1 
      NBLOCK = ( N - 1 ) / BLOCK_SIZE
      DO IBLOCK = 0, NBLOCK
        ROW_PROC = mod( IBLOCK, NPROW ) 
        IF ( MYROW.EQ.ROW_PROC ) THEN
          COL_PROC = mod( IBLOCK, NPCOL )
          IF ( MYCOL.EQ.COL_PROC ) THEN
            ILOC = ( IBLOCK / NPROW ) * BLOCK_SIZE
            JLOC = ( IBLOCK / NPCOL ) * BLOCK_SIZE
            I   =   ILOC + JLOC *  LOCAL_M + 1
            IMX = min(ILOC+BLOCK_SIZE,LOCAL_M)
     &            + (min(JLOC+BLOCK_SIZE,LOCAL_N)-1)*LOCAL_M
     &            + 1
            K=1
            DO WHILE ( I .LT. IMX )
             IF (SYM.NE.1) THEN
              ABSPIVOT =  abs(A(I))
             ELSE
              ABSPIVOT =  abs(A(I)*A(I))
             ENDIF
             CALL EMUMPS_UPDATE_MINMAX_PIVOT             ( ABSPIVOT,    
     +          DKEEP, KEEP, .FALSE.)
              K = K + 1
              I = I + DI
            END DO
          END IF
        END IF
      END DO
      RETURN
      END SUBROUTINE EMUMPS_PAR_ROOT_MINMAX_PIV_UPD
