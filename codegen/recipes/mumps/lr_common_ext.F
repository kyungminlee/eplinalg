C
C  Extended-precision MUMPS_LRMIXED_* helper stubs (kind16: Q / X).
C
C  MUMPS 5.9.0 added mixed-precision low-rank (LRMIXED) storage: BLR blocks
C  are optionally compressed into reduced IEEE byte formats (fp16/fp24/...).
C  lr_common.F hardcodes the CAST/PLUS/GEMM helpers for the genuine s/d/c/z
C  arithmetics; the migrated ?lr_core.F dispatchers for the extended prefixes
C  call renamed variants (CAST_QD, QGEMM, PLUS_...) that have no genuine
C  definition.  The feature is gated at every call site by
C  ``KEEP(448).GE.1 .AND. KEEP(450).NE.0`` and is never enabled for extended
C  precision (the byte-packed formats have no meaning for __float128 /
C  long double / double-double).  These stubs resolve the references and abort
C  loudly if the feature is ever forced on for an extended arithmetic.
C
C  Each stub inlines CALL MUMPS_ABORT() rather than sharing a helper, so no
C  symbol is defined by more than one arithmetic tree (all three trees'
C  precision archives link into a single mmsolve).
C
      SUBROUTINE MUMPS_LRMIXED_CAST_QD(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_QD: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_QD
      SUBROUTINE MUMPS_LRMIXED_CAST_QS(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_QS: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_QS
      SUBROUTINE MUMPS_LRMIXED_CAST_Q16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_Q16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_Q16
      SUBROUTINE MUMPS_LRMIXED_CAST_Q24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_Q24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_Q24
      SUBROUTINE MUMPS_LRMIXED_CAST_16Q(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16Q: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16Q
      SUBROUTINE MUMPS_LRMIXED_CAST_24Q(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24Q: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24Q
      SUBROUTINE MUMPS_LRMIXED_PLUS_QD(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_QD: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_QD
      SUBROUTINE MUMPS_LRMIXED_PLUS_QS(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_QS: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_QS
      SUBROUTINE MUMPS_LRMIXED_QGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_QGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_QGEMM
      SUBROUTINE MUMPS_LRMIXED_CAST_XC(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_XC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_XC
      SUBROUTINE MUMPS_LRMIXED_CAST_XZ(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_XZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_XZ
      SUBROUTINE MUMPS_LRMIXED_CAST_X16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_X16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_X16
      SUBROUTINE MUMPS_LRMIXED_CAST_X24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_X24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_X24
      SUBROUTINE MUMPS_LRMIXED_CAST_16X(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16X: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16X
      SUBROUTINE MUMPS_LRMIXED_CAST_24X(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24X: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24X
      SUBROUTINE MUMPS_LRMIXED_PLUS_XC(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_XC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_XC
      SUBROUTINE MUMPS_LRMIXED_PLUS_XZ(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_XZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_XZ
      SUBROUTINE MUMPS_LRMIXED_XGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_XGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_XGEMM
