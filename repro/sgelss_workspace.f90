program sgelss_workspace
   ! Reproducer for LAPACK 3.12.1 sgelss.f workspace-query
   ! inconsistency.
   !
   ! sgelss.f line 324 (path 2a, M < N with N >= MNTHR) computes the
   ! optimal workspace for its inner SGELQF call by:
   !
   !     MAXWRK = M + M * ILAENV(1, 'SGELQF', ' ', M, N, -1, -1)
   !
   ! while dgelss.f, in the analogous spot (lines 308-310 and 328),
   ! queries DGELQF directly with LWORK=-1 and uses the returned size:
   !
   !     CALL DGELQF(M, N, A, LDA, DUM, DUM, -1, INFO)
   !     LWORK_DGELQF = INT(DUM(1))
   !     ...
   !     MAXWRK = M + LWORK_DGELQF
   !
   ! Both forms happen to agree for the reference Netlib SGELQF (its
   ! optimal LWORK equals M*NB exactly), so there is no crash on
   ! reference builds — only a textual divergence vs the D-half.
   ! The mismatch becomes consequential if a vendor SGELQF reports
   ! an optimal LWORK that includes any constant beyond M*NB
   ! (workspace for the T factor, lookahead block, etc.).
   !
   ! This program demonstrates the inconsistency: it queries SGELSS
   ! and SGELQF independently, prints what each reports, and shows
   ! that the corresponding D-half SGELSS/DGELQF query pair is
   ! internally consistent in the same way.

   implicit none

   ! Pick M, N, NRHS that route SGELSS into path 2a (M < N,
   ! N >= MNTHR(SGELSS, M, N, NRHS)). For Netlib's ILAENV defaults
   ! MNTHR ~= 1.6*M, so N = 8*M is comfortably above the threshold.
   integer, parameter :: M = 64, N = 512, NRHS = 1
   integer, parameter :: LDA = M, LDB = N

   real,             allocatable :: SA(:,:), SB(:,:), SS(:)
   double precision, allocatable :: DA(:,:), DB(:,:), DS(:)

   real             :: SDUM(1)
   double precision :: DDUM(1)
   integer          :: info, lwork_sgelss, lwork_dgelss
   integer          :: lwork_sgelqf, lwork_dgelqf, rcond_int
   real             :: srcond
   double precision :: drcond

   external SGELSS, SGELQF, DGELSS, DGELQF
   integer  ILAENV
   external ILAENV

   allocate(SA(LDA, N), SB(LDB, NRHS), SS(min(M,N)))
   allocate(DA(LDA, N), DB(LDB, NRHS), DS(min(M,N)))
   SA = 0.0
   DA = 0.0d0
   SB = 0.0
   DB = 0.0d0
   srcond = -1.0
   drcond = -1.0d0

   ! Query SGELSS for its optimal workspace.
   call SGELSS(M, N, NRHS, SA, LDA, SB, LDB, SS, srcond, rcond_int, &
               SDUM, -1, info)
   if (info /= 0) then
      write(*,*) 'SGELSS query returned INFO=', info
      stop 1
   end if
   lwork_sgelss = int(SDUM(1))

   ! Query DGELSS for its optimal workspace.
   call DGELSS(M, N, NRHS, DA, LDA, DB, LDB, DS, drcond, rcond_int, &
               DDUM, -1, info)
   if (info /= 0) then
      write(*,*) 'DGELSS query returned INFO=', info
      stop 1
   end if
   lwork_dgelss = int(DDUM(1))

   ! Query SGELQF / DGELQF directly with the same M, N. This is the
   ! workspace the inner call inside the SGELSS / DGELSS path 2a
   ! actually needs.
   call SGELQF(M, N, SA, LDA, SDUM, SDUM, -1, info)
   if (info /= 0) then
      write(*,*) 'SGELQF query returned INFO=', info
      stop 1
   end if
   lwork_sgelqf = int(SDUM(1))

   call DGELQF(M, N, DA, LDA, DDUM, DDUM, -1, info)
   if (info /= 0) then
      write(*,*) 'DGELQF query returned INFO=', info
      stop 1
   end if
   lwork_dgelqf = int(DDUM(1))

   write(*,'(A)') 'sgelss workspace reproducer:'
   write(*,'(A,I0,A,I0,A,I0)') '  M=', M, '  N=', N, '  NRHS=', NRHS

   write(*,*)
   write(*,'(A)') '=== sgelss.f path 2a SGELQF term ==='
   write(*,'(A,I0,A,I0)') &
        '  formula (line 324):  M + M*ILAENV(1,SGELQF) = ', &
        M, ' + ', M * ILAENV(1, 'SGELQF', ' ', M, N, -1, -1)
   write(*,'(A,I0,A,I0)') &
        '                                    total =  ', &
        M + M * ILAENV(1, 'SGELQF', ' ', M, N, -1, -1)
   write(*,'(A,I0,A,I0)') &
        '  direct  (would-be):  M + LWORK_SGELQF      = ', &
        M, ' + ', lwork_sgelqf
   write(*,'(A,I0)') &
        '                                    total =  ', &
        M + lwork_sgelqf
   write(*,*)
   write(*,'(A)') '=== dgelss.f path 2a DGELQF term (already correct) ==='
   write(*,'(A,I0,A,I0)') &
        '  formula (line 328):  M + LWORK_DGELQF      = ', &
        M, ' + ', lwork_dgelqf
   write(*,'(A,I0)') &
        '                                    total =  ', &
        M + lwork_dgelqf
   write(*,*)
   write(*,'(A,I0)') '  SGELSS reported optimal LWORK = ', lwork_sgelss
   write(*,'(A,I0)') '  DGELSS reported optimal LWORK = ', lwork_dgelss
   write(*,*)
   write(*,'(A)') '=== Verdict ==='
   if (M * ILAENV(1, 'SGELQF', ' ', M, N, -1, -1) == lwork_sgelqf) then
      write(*,'(A)') '  Reference SGELQF: ILAENV formula and direct'
      write(*,'(A)') '  query coincide (both = M*NB). Bug is LATENT'
      write(*,'(A)') '  in the reference build — no visible crash.'
   else
      write(*,'(A,I0,A,I0,A)') &
           '  Reference SGELQF: ILAENV formula (', &
           M * ILAENV(1, 'SGELQF', ' ', M, N, -1, -1), &
           ') != direct (', lwork_sgelqf, ')  ⚠'
   end if
   write(*,*)
   write(*,'(A)') '  The S-half STILL uses the indirect ILAENV form'
   write(*,'(A)') '  (sgelss.f:324) instead of the direct query'
   write(*,'(A)') '  pattern that the D-half adopted (dgelss.f:307-310,'
   write(*,'(A)') '  328). Vendor or future-reference SGELQF whose'
   write(*,'(A)') '  optimal LWORK exceeds M*NB would silently'
   write(*,'(A)') '  underestimate, leading to insufficient workspace'
   write(*,'(A)') '  on the SGELSS path 2a branch.'

   deallocate(SA, SB, SS, DA, DB, DS)
end program sgelss_workspace
