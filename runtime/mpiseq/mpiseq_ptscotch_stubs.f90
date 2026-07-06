! PT-Scotch dgraph stub forwarders for libmpiseq.
!
! When PT-Scotch is staged (real-MPI columns), the migrated mumps
! archives are compiled with -Dptscotch and their *ana_aux_par.F
! members reference the distributed-graph Fortran API
! (SCOTCHFDGRAPH*_MUMPS / SCOTCHFSTRATDGRAPHORDER_MUMPS) plus the C
! entry SCOTCH_dgraphInit_mumps via the MUMPS_DGRAPHINIT wrapper in
! mumps_scotch.c. Those symbols live in the real ptscotch archive,
! which cannot appear on a libmpiseq link line — PT-Scotch requires a
! real MPI underneath. This file ships aborting stand-ins so the
! ``_seq`` link of those same archives resolves; a single-rank
! libmpiseq executable can never reach distributed ordering (MUMPS
! forces sequential analysis on one process), so like the rest of the
! libmpiseq surface each stub prints "should not be called" and stops.
! Sequential Scotch (ICNTL(7)=3) is untouched — its symbols resolve
! from the ordinary scotch archive.
!
! The ``_MUMPS`` name suffix matches the vendored Scotch build
! (SCOTCH_NAME_SUFFIX=_mumps), which is the flavor the migrated mumps
! archives are compiled against. The stub list is exactly the surface
! *ana_aux_par.F references without MUMPS_SCOTCHIMPORTOMPTHREADS
! (SCOTCHFCONTEXT* are guarded behind that define upstream).
! Argument lists mirror libscotch's library_dgraph*_f.c so the arity
! survives cross-unit interface checks; no stub ever touches them.

subroutine SCOTCHFDGRAPHBUILD_MUMPS(grafdat, baseval, vertlocnbr, &
    vertlocmax, vertloctab, vendloctab, veloloctab, vlblloctab, &
    edgelocnbr, edgelocsiz, edgeloctab, edgegsttab, edloloctab, ierr)
  implicit none
  double precision :: grafdat(*)
  integer :: baseval, vertlocnbr, vertlocmax
  integer :: vertloctab(*), vendloctab(*), veloloctab(*), vlblloctab(*)
  integer :: edgelocnbr, edgelocsiz
  integer :: edgeloctab(*), edgegsttab(*), edloloctab(*)
  integer :: ierr
  write(*,*) 'Error. SCOTCHFDGRAPHBUILD should not be called.'
  stop
end subroutine SCOTCHFDGRAPHBUILD_MUMPS

subroutine SCOTCHFDGRAPHEXIT_MUMPS(grafdat)
  implicit none
  double precision :: grafdat(*)
  write(*,*) 'Error. SCOTCHFDGRAPHEXIT should not be called.'
  stop
end subroutine SCOTCHFDGRAPHEXIT_MUMPS

subroutine SCOTCHFDGRAPHORDERINIT_MUMPS(grafdat, ordedat, ierr)
  implicit none
  double precision :: grafdat(*), ordedat(*)
  integer :: ierr
  write(*,*) 'Error. SCOTCHFDGRAPHORDERINIT should not be called.'
  stop
end subroutine SCOTCHFDGRAPHORDERINIT_MUMPS

subroutine SCOTCHFDGRAPHORDEREXIT_MUMPS(grafdat, ordedat)
  implicit none
  double precision :: grafdat(*), ordedat(*)
  write(*,*) 'Error. SCOTCHFDGRAPHORDEREXIT should not be called.'
  stop
end subroutine SCOTCHFDGRAPHORDEREXIT_MUMPS

subroutine SCOTCHFDGRAPHORDERCOMPUTE_MUMPS(grafdat, ordedat, stradat, &
    ierr)
  implicit none
  double precision :: grafdat(*), ordedat(*), stradat(*)
  integer :: ierr
  write(*,*) 'Error. SCOTCHFDGRAPHORDERCOMPUTE should not be called.'
  stop
end subroutine SCOTCHFDGRAPHORDERCOMPUTE_MUMPS

subroutine SCOTCHFDGRAPHORDERGATHER_MUMPS(grafdat, dorddat, corddat, &
    ierr)
  implicit none
  double precision :: grafdat(*), dorddat(*), corddat(*)
  integer :: ierr
  write(*,*) 'Error. SCOTCHFDGRAPHORDERGATHER should not be called.'
  stop
end subroutine SCOTCHFDGRAPHORDERGATHER_MUMPS

subroutine SCOTCHFDGRAPHCORDERINIT_MUMPS(grafdat, corddat, permtab, &
    peritab, cblknbr, rangtab, treetab, ierr)
  implicit none
  double precision :: grafdat(*), corddat(*)
  integer :: permtab(*), peritab(*), rangtab(*), treetab(*)
  integer :: cblknbr, ierr
  write(*,*) 'Error. SCOTCHFDGRAPHCORDERINIT should not be called.'
  stop
end subroutine SCOTCHFDGRAPHCORDERINIT_MUMPS

subroutine SCOTCHFDGRAPHCORDEREXIT_MUMPS(grafdat, corddat)
  implicit none
  double precision :: grafdat(*), corddat(*)
  write(*,*) 'Error. SCOTCHFDGRAPHCORDEREXIT should not be called.'
  stop
end subroutine SCOTCHFDGRAPHCORDEREXIT_MUMPS

subroutine SCOTCHFSTRATDGRAPHORDER_MUMPS(stradat, string, ierr)
  implicit none
  double precision :: stradat(*)
  character(len=*) :: string
  integer :: ierr
  write(*,*) 'Error. SCOTCHFSTRATDGRAPHORDER should not be called.'
  stop
end subroutine SCOTCHFSTRATDGRAPHORDER_MUMPS

! C-side twin: mumps_scotch.c's MUMPS_DGRAPHINIT calls
! SCOTCH_dgraphInit(graf, MPI_Comm_f2c(*comm)) — with the vendored
! rename that is ``int SCOTCH_dgraphInit_mumps(SCOTCH_Dgraph *,
! MPI_Comm)``, both arguments by value (Intel MPI's MPI_Comm is int).
integer(c_int) function SCOTCH_dgraphInit_mumps(grafptr, comm) &
    bind(c, name='SCOTCH_dgraphInit_mumps')
  use iso_c_binding, only: c_int, c_ptr
  implicit none
  type(c_ptr), value :: grafptr
  integer(c_int), value :: comm
  write(*,*) 'Error. SCOTCH_dgraphInit should not be called.'
  SCOTCH_dgraphInit_mumps = 1
  stop
end function SCOTCH_dgraphInit_mumps
