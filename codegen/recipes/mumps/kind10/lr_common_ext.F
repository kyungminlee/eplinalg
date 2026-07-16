C
C  Extended-precision MUMPS_LRMIXED_* helper stubs (kind10: E / Y).
C
C  MUMPS 5.9.0 added mixed-precision low-rank (LRMIXED) storage: BLR blocks
C  are optionally compressed into reduced IEEE byte formats (fp16/fp24/...).
C  lr_common.F hardcodes the CAST/PLUS/GEMM helpers for the genuine s/d/c/z
C  arithmetics; the migrated ?lr_core.F dispatchers for the extended prefixes
C  call renamed variants (CAST_QD, EGEMM, PLUS_...) that have no genuine
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
      SUBROUTINE MUMPS_LRMIXED_CAST_ED(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_ED: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_ED
      SUBROUTINE MUMPS_LRMIXED_CAST_ES(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_ES: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_ES
      SUBROUTINE MUMPS_LRMIXED_CAST_E16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_E16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_E16
      SUBROUTINE MUMPS_LRMIXED_CAST_E24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_E24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_E24
      SUBROUTINE MUMPS_LRMIXED_CAST_16E(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16E: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16E
      SUBROUTINE MUMPS_LRMIXED_CAST_24E(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24E: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24E
      SUBROUTINE MUMPS_LRMIXED_PLUS_ED(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_ED: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_ED
      SUBROUTINE MUMPS_LRMIXED_PLUS_ES(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_ES: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_ES
      SUBROUTINE MUMPS_LRMIXED_EGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_EGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_EGEMM
      SUBROUTINE MUMPS_LRMIXED_CAST_YC(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_YC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_YC
      SUBROUTINE MUMPS_LRMIXED_CAST_YZ(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_YZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_YZ
      SUBROUTINE MUMPS_LRMIXED_CAST_Y16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_Y16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_Y16
      SUBROUTINE MUMPS_LRMIXED_CAST_Y24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_Y24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_Y24
      SUBROUTINE MUMPS_LRMIXED_CAST_16Y(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16Y: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16Y
      SUBROUTINE MUMPS_LRMIXED_CAST_24Y(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24Y: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24Y
      SUBROUTINE MUMPS_LRMIXED_PLUS_YC(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_YC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_YC
      SUBROUTINE MUMPS_LRMIXED_PLUS_YZ(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_YZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_YZ
      SUBROUTINE MUMPS_LRMIXED_YGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_YGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_YGEMM
