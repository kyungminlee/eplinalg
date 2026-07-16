C
C  Hand-written extended-precision reallocator module.
C
C  Bridges callers migrated to extended-precision arithmetic (Q / X
C  for kind16; extendable to E / Y / DD / ZZ on platforms that
C  support them) to the reallocator contract provided by the
C  upstream MUMPS_MEMORY_MOD.
C
C  The original MUMPS_MEMORY_MOD bundles MUMPS_S/D/C/ZREALLOC under
C  a single MUMPS_REALLOC generic interface. When the migrator
C  collapses S and D to Q (and C and Z to X), those generic bundles
C  produce duplicate MODULE PROCEDURE entries. We dodge that by
C  keeping the upstream module unchanged (copy_files) and providing
C  a separate module that re-exports the original procedures and
C  contributes the new extended-precision ones.
C
C  The migrator rewrites ``USE MUMPS_MEMORY_MOD`` → ``USE
C  MUMPS_MEMORY_MOD_EP`` in migrated callers via the recipe's
C  ``module_renames`` field; callers that are copied verbatim
C  continue to use MUMPS_MEMORY_MOD directly.
C
      MODULE MUMPS_MEMORY_MOD_EP
C     Import only the individual reallocator subroutines from the
C     upstream module; do NOT import its ``MUMPS_REALLOC`` /
C     ``MUMPS_DEALLOC`` generic interfaces — we redefine them below
C     to cover the extended-precision variants alongside the originals.
      USE MUMPS_MEMORY_MOD, ONLY:
     &     MUMPS_IREALLOC,  MUMPS_IREALLOC8,
     &     MUMPS_I8REALLOC, MUMPS_I8REALLOC8,
     &     MUMPS_SREALLOC,  MUMPS_DREALLOC,
     &     MUMPS_CREALLOC,  MUMPS_ZREALLOC,
     &     MUMPS_IDEALLOC,  MUMPS_I8DEALLOC,
     &     MUMPS_MEMORY_SET_DATA_SIZES
      IMPLICIT NONE
      INTERFACE MUMPS_REALLOC
      MODULE PROCEDURE MUMPS_IREALLOC
      MODULE PROCEDURE MUMPS_SREALLOC, MUMPS_DREALLOC
      MODULE PROCEDURE MUMPS_CREALLOC, MUMPS_ZREALLOC
      MODULE PROCEDURE MUMPS_QREALLOC, MUMPS_XREALLOC
      END INTERFACE
      INTERFACE MUMPS_DEALLOC
      MODULE PROCEDURE MUMPS_IDEALLOC
      END INTERFACE
      CONTAINS
      SUBROUTINE MUMPS_QREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      REAL(KIND=16), POINTER       :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      REAL(KIND=16), POINTER       :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      INTEGER(8)                   :: ELSIZE
      CHARACTER(len=60)            :: FMTA, FMTD
      ELSIZE = int(storage_size(1.0_16),8)/8_8
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*ELSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+int(MINSIZE,8)*ELSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_QREALLOC
      SUBROUTINE MUMPS_XREALLOC(ARRAY, MINSIZE, INFO, LP, FORCE, COPY,
     &     STRING, MEMCNT, ERRCODE)
      COMPLEX(KIND=16), POINTER    :: ARRAY(:)
      INTEGER                      :: INFO(:)
      INTEGER                      :: MINSIZE, LP
      LOGICAL, OPTIONAL            :: FORCE
      LOGICAL, OPTIONAL            :: COPY
      CHARACTER, OPTIONAL          :: STRING*(*)
      INTEGER, OPTIONAL            :: ERRCODE
      INTEGER(8), OPTIONAL         :: MEMCNT
      LOGICAL                      :: ICOPY, IFORCE
      COMPLEX(KIND=16), POINTER    :: TEMP(:)
      INTEGER                      :: I, IERR, ERRTPL(2)
      INTEGER(8)                   :: ELSIZE
      CHARACTER(len=60)            :: FMTA, FMTD
      ELSIZE = int(storage_size((1.0_16,0.0_16)),8)/8_8
      IF(present(COPY)) THEN
         ICOPY = COPY
      ELSE
         ICOPY = .FALSE.
      END IF
      IF (present(FORCE)) THEN
         IFORCE = FORCE
      ELSE
         IFORCE = .FALSE.
      END IF
      IF (present(STRING)) THEN
         FMTA = "Allocation failed inside realloc: "//STRING
         FMTD = "Deallocation failed inside realloc: "//STRING
      ELSE
         FMTA = "Allocation failed inside realloc: "
         FMTD = "Deallocation failed inside realloc: "
      END IF
      IF (present(ERRCODE)) THEN
         ERRTPL(1) = ERRCODE
         ERRTPL(2) = MINSIZE
      ELSE
         ERRTPL(1) = -13
         ERRTPL(2) = MINSIZE
      END IF
      IF(ICOPY) THEN
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               allocate(TEMP(MINSIZE), STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTA)
                  INFO(1:2) = ERRTPL
                  RETURN
               ELSE
                  IF(present(MEMCNT))MEMCNT = MEMCNT+
     &                 int(MINSIZE,8)*ELSIZE
               END IF
               DO I=1, min(size(ARRAY), MINSIZE)
                  TEMP(I) = ARRAY(I)
               END DO
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
               NULLIFY(ARRAY)
               ARRAY => TEMP
               NULLIFY(TEMP)
            END IF
         ELSE
            WRITE(LP,
     &      '("Input array is not associated. nothing to copy here")')
            RETURN
         END IF
      ELSE
         IF(associated(ARRAY)) THEN
            IF ((size(ARRAY) .LT. MINSIZE) .OR.
     &           ((size(ARRAY).NE.MINSIZE) .AND. IFORCE)) THEN
               IF(present(MEMCNT))MEMCNT = MEMCNT-
     &              int(size(ARRAY),8)*ELSIZE
               deallocate(ARRAY, STAT=IERR)
               IF(IERR .LT. 0) THEN
                  WRITE(LP,FMTD)
                  INFO(1:2) = ERRTPL
                  RETURN
               END IF
            ELSE
               RETURN
            END IF
         END IF
         allocate(ARRAY(MINSIZE), STAT=IERR)
         IF(IERR .LT. 0) THEN
            WRITE(LP,FMTA)
            INFO(1:2) = ERRTPL
            RETURN
         ELSE
            IF(present(MEMCNT)) MEMCNT = MEMCNT+int(MINSIZE,8)*ELSIZE
         END IF
      END IF
      RETURN
      END SUBROUTINE MUMPS_XREALLOC
      END MODULE MUMPS_MEMORY_MOD_EP
