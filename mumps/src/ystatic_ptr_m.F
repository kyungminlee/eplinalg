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
      MODULE YMUMPS_STATIC_PTR_M
      PUBLIC :: ZMUMPS_TMP_PTR, YMUMPS_GET_TMP_PTR
      COMPLEX(KIND=10), DIMENSION(:), POINTER, SAVE :: ZMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE YMUMPS_SET_STATIC_PTR(ARRAY)
      COMPLEX(KIND=10), DIMENSION(:), TARGET :: ARRAY
      ZMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE YMUMPS_SET_STATIC_PTR
      SUBROUTINE YMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_NOF2003)
      COMPLEX(KIND=10), DIMENSION(:), POINTER :: PTR
#else
      COMPLEX(KIND=10), DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#endif
      PTR => ZMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE YMUMPS_GET_TMP_PTR
      END MODULE YMUMPS_STATIC_PTR_M
