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
      MODULE EMUMPS_STATIC_PTR_M
      PUBLIC :: DMUMPS_TMP_PTR, EMUMPS_GET_TMP_PTR
      REAL(KIND=10), DIMENSION(:), POINTER, SAVE :: DMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE EMUMPS_SET_STATIC_PTR(ARRAY)
      REAL(KIND=10), DIMENSION(:), TARGET :: ARRAY
      DMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE EMUMPS_SET_STATIC_PTR
      SUBROUTINE EMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_NOF2003)
      REAL(KIND=10), DIMENSION(:), POINTER :: PTR
#else
      REAL(KIND=10), DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#endif
      PTR => DMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE EMUMPS_GET_TMP_PTR
      END MODULE EMUMPS_STATIC_PTR_M
