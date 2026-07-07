C
C  Extended-precision MUMPS_LRMIXED_* helper stubs (multifloats: M / W).
C
C  MUMPS 5.9.0 added mixed-precision low-rank (LRMIXED) storage: BLR blocks
C  are optionally compressed into reduced IEEE byte formats (fp16/fp24/...).
C  lr_common.F hardcodes the CAST/PLUS/GEMM helpers for the genuine s/d/c/z
C  arithmetics; the migrated ?lr_core.F dispatchers for the extended prefixes
C  call renamed variants (CAST_QD, MGEMM, PLUS_...) that have no genuine
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
      SUBROUTINE MUMPS_LRMIXED_CAST_MD(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_MD: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_MD
      SUBROUTINE MUMPS_LRMIXED_CAST_MS(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_MS: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_MS
      SUBROUTINE MUMPS_LRMIXED_CAST_M16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_M16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_M16
      SUBROUTINE MUMPS_LRMIXED_CAST_M24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_M24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_M24
      SUBROUTINE MUMPS_LRMIXED_CAST_16M(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16M: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16M
      SUBROUTINE MUMPS_LRMIXED_CAST_24M(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24M: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24M
      SUBROUTINE MUMPS_LRMIXED_PLUS_MD(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_MD: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_MD
      SUBROUTINE MUMPS_LRMIXED_PLUS_MS(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_MS: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_MS
      SUBROUTINE MUMPS_LRMIXED_MGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_MGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_MGEMM
      SUBROUTINE MUMPS_LRMIXED_CAST_WC(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_WC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_WC
      SUBROUTINE MUMPS_LRMIXED_CAST_WZ(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_WZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_WZ
      SUBROUTINE MUMPS_LRMIXED_CAST_W16(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_W16: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_W16
      SUBROUTINE MUMPS_LRMIXED_CAST_W24(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_W24: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_W24
      SUBROUTINE MUMPS_LRMIXED_CAST_16W(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_16W: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_16W
      SUBROUTINE MUMPS_LRMIXED_CAST_24W(A, B, NROWS, NCOLS, LDA, LDB)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDA, LDB
      INTEGER :: A(*), B(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_CAST_24W: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_CAST_24W
      SUBROUTINE MUMPS_LRMIXED_PLUS_WC(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_WC: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_WC
      SUBROUTINE MUMPS_LRMIXED_PLUS_WZ(X, Y, NROWS, NCOLS, LDX, LDY)
      IMPLICIT NONE
      INTEGER :: NROWS, NCOLS, LDX, LDY
      INTEGER :: X(*), Y(*)
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_PLUS_WZ: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_PLUS_WZ
      SUBROUTINE MUMPS_LRMIXED_WGEMM(TRANSA, TRANSB, M, N, K, ALPHA,
     &               A, LDA, B, LDB, BETA, C, LDC)
      IMPLICIT NONE
      CHARACTER(len=1) :: TRANSA, TRANSB
      INTEGER :: M, N, K, LDA, LDB, LDC
      INTEGER :: A(*), B(*), C(*), ALPHA, BETA
      WRITE(*,*) 'FATAL: MUMPS_LRMIXED_WGEMM: mixed-precision',
     &  ' low-rank (KEEP(450)/=0) unsupported for extended precision.'
      CALL MUMPS_ABORT()
      END SUBROUTINE MUMPS_LRMIXED_WGEMM
