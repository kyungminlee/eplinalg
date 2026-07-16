"""MUMPS libseq ``mpi.f`` patch — a recipe-specific staging hack.

Extracted verbatim from ``__main__.py`` (Cluster 3) as part of the migrator
file-restructuring refactor. Behaviour is unchanged; ``cmd_stage`` imports
``_patch_libseq_mpi_f`` from here.
"""
from pathlib import Path


def _patch_libseq_mpi_f(path: Path) -> None:
    """Extend libseq's ``MUMPS_COPY`` with MPI_REAL16 / MPI_COMPLEX32
    cases so reductions on REAL(KIND=16) / COMPLEX(KIND=16) buffers
    dispatch correctly under our libmpiseq variant. Upstream's
    ``MUMPS_COPY`` only knows the standard MPI datatypes; the migrated
    qxmumps archive passes MPI_REAL16 (Intel MPI = 1275072555) for
    kind16 reductions.

    Patches the staged copy at ``_mpiseq_src/mpi.f``; upstream's
    ``extern/MUMPS_5.9.0/libseq/mpi.f`` stays read-only. BLACS /
    ScaLAPACK forwarders inside the same file are deliberately KEPT
    — libmpiseq stands in for those archives in the ``_seq`` test
    link, and the real BLACS / ScaLAPACK archives aren't linked there
    so there's no duplicate-symbol collision.
    """
    src = path.read_text()

    # Extend MUMPS_COPY's dispatch with MPI_REAL16 / MPI_COMPLEX32
    # cases, and append matching MUMPS_COPY_* helpers. Anchor the
    # insertion before the existing ``ELSE\n IERR=1`` fallthrough.
    fallthrough = '      ELSE\n        IERR=1\n        RETURN\n      END IF'
    extra_dispatch = (
        # kind16: REAL(16) / COMPLEX(16)
        '      ELSE IF ( DATATYPE .EQ. MPI_REAL16 ) THEN\n'
        '      CALL MUMPS_COPY_REAL16( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_COMPLEX32 ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX32( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # kind10: 80-bit extended real / complex map to MPI's long
        # double tokens (no MPI_REAL10 in standard MPI).
        '      ELSE IF ( DATATYPE .EQ. MPI_LONG_DOUBLE ) THEN\n'
        '      CALL MUMPS_COPY_REAL10( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. MPI_C_LONG_DOUBLE_COMPLEX ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX20( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        # multifloats: runtime/mpiseq/mpiseq_c_stubs.c encodes derived-type
        # sentinels as 0x10000000 | total_bytes. float64x2 → 16-byte
        # element (sentinel 268435472); complex64x2 → 32-byte
        # (268435488). MPI_Type_c2f / MPI_Op_c2f in Intel mpi.h are
        # passthrough casts, so the Fortran handle is the same value.
        '      ELSE IF ( DATATYPE .EQ. 268435472 ) THEN\n'
        '      CALL MUMPS_COPY_FLOAT64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
        '      ELSE IF ( DATATYPE .EQ. 268435488 ) THEN\n'
        '      CALL MUMPS_COPY_COMPLEX64X2( SENDBUF, RECVBUF, CNT, SS, RS )\n'
    )
    if 'MPI_REAL16' not in src and fallthrough in src:
        src = src.replace(fallthrough, extra_dispatch + fallthrough, 1)

    extra_helpers = """
      SUBROUTINE MUMPS_COPY_REAL16( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=16) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_REAL16
      SUBROUTINE MUMPS_COPY_COMPLEX32( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      COMPLEX(KIND=16) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX32
      SUBROUTINE MUMPS_COPY_REAL10( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=10) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_REAL10
      SUBROUTINE MUMPS_COPY_COMPLEX20( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      COMPLEX(KIND=10) S(N),R(N)
      INTEGER I
      DO I = 1, N
        R(I+RS) = S(I+SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX20
      SUBROUTINE MUMPS_COPY_FLOAT64X2( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=8) S(2*N),R(2*N)
      INTEGER I
      DO I = 1, 2*N
        R(I+2*RS) = S(I+2*SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_FLOAT64X2
      SUBROUTINE MUMPS_COPY_COMPLEX64X2( S, R, N, SS, RS )
      IMPLICIT NONE
      INTEGER N, SS, RS
      REAL(KIND=8) S(4*N),R(4*N)
      INTEGER I
      DO I = 1, 4*N
        R(I+4*RS) = S(I+4*SS)
      END DO
      RETURN
      END SUBROUTINE MUMPS_COPY_COMPLEX64X2
"""
    if 'SUBROUTINE MUMPS_COPY_REAL16' not in src:
        src = src.rstrip() + '\n' + extra_helpers

    # Base libseq mpi.f stubs pchk2mat (MUMPS calls it) but NOT pchk1mat /
    # globchk. Our migrated *typed* ScaLAPACK archives (ey/qx/mw) are full
    # netlib ports that reference all three, co-located in the reference
    # ScaLAPACK's single pchkxmat.f.o object. In a seq link, pulling that
    # object in for pchk1mat_/globchk_ also drags its pchk2mat_, colliding
    # with libseq's pchk2mat_. Stubbing pchk1mat / globchk here makes libseq
    # fully cover pchkxmat.f.o, so the reference object is never extracted and
    # there is no duplicate symbol (no --allow-multiple-definition). Like
    # pchk2mat these are descriptor checks reached only from inside pXgetrf/
    # pXpotrf, which MUMPS never calls at np=1 (the root is factored by
    # sequential LAPACK) — so the STOP body never executes.
    extra_pchk = """
      SUBROUTINE pchk1mat( MA, MAPOS0, NA, NAPOS0, IA, JA, DESCA,
     &                     DESCAPOS0, NEXTRA, EX, EXPOS, INFO )
      IMPLICIT NONE
      INTEGER            DESCAPOS0, IA, INFO, JA, MA, MAPOS0, NA,
     &                   NAPOS0, NEXTRA
      INTEGER            DESCA( * ), EX( NEXTRA ), EXPOS( NEXTRA )
        WRITE(*,*) 'Error. PCHK1MAT should not be called.'
        STOP
      RETURN
      END SUBROUTINE pchk1mat
      SUBROUTINE globchk( ICTXT, N, X, LDX, IWORK, INFO )
      IMPLICIT NONE
      INTEGER            ICTXT, INFO, LDX, N
      INTEGER            IWORK( * ), X( LDX, * )
        WRITE(*,*) 'Error. GLOBCHK should not be called.'
        STOP
      RETURN
      END SUBROUTINE globchk
"""
    if 'SUBROUTINE pchk1mat' not in src:
        src = src.rstrip() + '\n' + extra_pchk

    path.write_text(src)
