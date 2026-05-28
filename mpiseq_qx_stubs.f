C  Per-precision ScaLAPACK stub forwarders for libmpiseq.
C
C  Upstream's external/MUMPS_5.8.2/libseq/mpi.f provides single-process
C  "should not be called" stubs only for the D-prefixed ScaLAPACK
C  symbols (PDGETRF / PDGETRS / PDPOTRF / PDPOTRS / PDTRTRS). The
C  migrated qmumps / xmumps (kind16) and emumps / ymumps (kind10)
C  archives request the Q/X/E/Y-prefixed equivalents — those resolve
C  through ${LIB_PREFIX}scalapack at the standard mpiexec path, but
C  for a fully-sequential libmpiseq link those symbols have to live
C  somewhere too. This file ships them inside the migrator (the
C  external/ tree is read-only) and is appended to the mpiseq target
C  in cmake/CMakeLists.txt.
C
C  Each stub mirrors upstream's pattern exactly: print "should not be
C  called", STOP. Single-rank operation never reaches them; they exist
C  only so the link resolves.
C
C  Multifloats prefixes (M, W) live in mpiseq_mw_stubs.f90 — kept in a
C  separate free-form file because the stubs require TYPE(real64x2) /
C  TYPE(cmplx64x2) from the multifloats Fortran module.
C
C***********************************************************************
      SUBROUTINE PQGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL*16            A( * )
        WRITE(*,*) 'Error. PQGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PQGETRF
C***********************************************************************
      SUBROUTINE PQGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL*16            A( * ), B( * )
        WRITE(*,*) 'Error. PQGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PQGETRS
C***********************************************************************
      SUBROUTINE PQPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL*16            A( * )
        WRITE(*,*) 'Error. PQPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PQPOTRF
C***********************************************************************
      SUBROUTINE PQPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*16            A( * ), B( * )
        WRITE(*,*) 'Error. PQPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PQPOTRS
C***********************************************************************
      SUBROUTINE PQTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*16            A( * ), B( * )
        WRITE(*,*) 'Error. PQTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PQTRTRS
C***********************************************************************
      SUBROUTINE PXGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX*32         A( * )
        WRITE(*,*) 'Error. PXGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PXGETRF
C***********************************************************************
      SUBROUTINE PXGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*32         A( * ), B( * )
        WRITE(*,*) 'Error. PXGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PXGETRS
C***********************************************************************
      SUBROUTINE PXPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX*32         A( * )
        WRITE(*,*) 'Error. PXPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PXPOTRF
C***********************************************************************
      SUBROUTINE PXPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*32         A( * ), B( * )
        WRITE(*,*) 'Error. PXPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PXPOTRS
C***********************************************************************
      SUBROUTINE PXTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*32         A( * ), B( * )
        WRITE(*,*) 'Error. PXTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PXTRTRS
C***********************************************************************
      SUBROUTINE PEGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      REAL*10            A( * )
        WRITE(*,*) 'Error. PEGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PEGETRF
C***********************************************************************
      SUBROUTINE PEGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      REAL*10            A( * ), B( * )
        WRITE(*,*) 'Error. PEGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PEGETRS
C***********************************************************************
      SUBROUTINE PEPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      REAL*10            A( * )
        WRITE(*,*) 'Error. PEPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PEPOTRF
C***********************************************************************
      SUBROUTINE PEPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*10            A( * ), B( * )
        WRITE(*,*) 'Error. PEPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PEPOTRS
C***********************************************************************
      SUBROUTINE PETRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      REAL*10            A( * ), B( * )
        WRITE(*,*) 'Error. PETRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PETRTRS
C***********************************************************************
      SUBROUTINE PYGETRF( M, N, A, IA, JA, DESCA, IPIV, INFO )
      IMPLICIT NONE
      INTEGER            IA, INFO, JA, M, N
      INTEGER            DESCA( * ), IPIV( * )
      COMPLEX*20         A( * )
        WRITE(*,*) 'Error. PYGETRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PYGETRF
C***********************************************************************
      SUBROUTINE PYGETRS( TRANS, N, NRHS, A, IA, JA, DESCA, IPIV,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          TRANS
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * ), IPIV( * )
      COMPLEX*20         A( * ), B( * )
        WRITE(*,*) 'Error. PYGETRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PYGETRS
C***********************************************************************
      SUBROUTINE PYPOTRF( UPLO, N, A, IA, JA, DESCA, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, INFO, JA, N
      INTEGER            DESCA( * )
      COMPLEX*20         A( * )
        WRITE(*,*) 'Error. PYPOTRF should not be called.'
        STOP
      RETURN
      END SUBROUTINE PYPOTRF
C***********************************************************************
      SUBROUTINE PYPOTRS( UPLO, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*20         A( * ), B( * )
        WRITE(*,*) 'Error. PYPOTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PYPOTRS
C***********************************************************************
      SUBROUTINE PYTRTRS( UPLO, TRANS, DIAG, N, NRHS, A, IA, JA, DESCA,
     &                    B, IB, JB, DESCB, INFO )
      IMPLICIT NONE
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            IA, IB, INFO, JA, JB, N, NRHS
      INTEGER            DESCA( * ), DESCB( * )
      COMPLEX*20         A( * ), B( * )
        WRITE(*,*) 'Error. PYTRTRS should not be called.'
        STOP
      RETURN
      END SUBROUTINE PYTRTRS
