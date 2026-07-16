! multifloats_mpi_f -- Fortran-side MPI handles for the real-precision
! multifloats datatype and the real + complex reduction operators
! registered by multifloats_mpi_init().
!
! Migrated multifloats Fortran code (e.g. mmumps) calls MPI directly
! with statements like ``CALL MPI_SEND(..., MPI_FLOAT64X2, ...)`` and
! ``CALL MPI_REDUCE(..., MPI_MM_SUM, ...)``. The handles are created in
! C++ at runtime, not by mpif.h, so Fortran needs an explicit module
! that exposes them. Each public name below is a default-INTEGER
! (matching MPI_Fint on Linux/x86_64) bound to the C symbol
! multifloats_mpi.cpp populates via MPI_Type_c2f / MPI_Op_c2f. The
! reduction-op handles are 0 until ``multifloats_mpi_init`` runs;
! trigger it explicitly via ``CALL multifloats_mpi_init`` after
! MPI_Init for consumers (like MUMPS) under real MPI. The two datatype
! handles (MPI_FLOAT64X2 / MPI_COMPLEX64X2) instead default to the
! libmpiseq derived-type sentinel, so a sequential (libmpiseq) MUMPS
! consumer works with neither MPI_Init nor multifloats_mpi_init() — see
! the note in multifloats_mpi.cpp. Real-MPI init overwrites them.
!
module multifloats_mpi_f
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   private

   public :: MPI_FLOAT64X2, MPI_COMPLEX64X2
   public :: MPI_MM_SUM, MPI_MM_AMX, MPI_MM_AMN
   public :: MPI_WW_SUM, MPI_WW_AMX, MPI_WW_AMN
   public :: multifloats_mpi_init

   ! Deprecated pre-v0.14 names (arithmetic-based DD/ZZ); use the
   ! family-prefixed MPI_MM_* / MPI_WW_* above. Bound to the C alias
   ! symbols of the same handles. Remove after v0.14.
   public :: MPI_DD_SUM, MPI_DD_AMX, MPI_DD_AMN
   public :: MPI_ZZ_SUM, MPI_ZZ_AMX, MPI_ZZ_AMN

   integer(c_int), bind(c, name='mf_mpi_float64x2_f'),  protected :: MPI_FLOAT64X2
   integer(c_int), bind(c, name='mf_mpi_complex64x2_f'),protected :: MPI_COMPLEX64X2
   integer(c_int), bind(c, name='mf_mpi_mm_sum_f'),     protected :: MPI_MM_SUM
   integer(c_int), bind(c, name='mf_mpi_mm_amx_f'),     protected :: MPI_MM_AMX
   integer(c_int), bind(c, name='mf_mpi_mm_amn_f'),     protected :: MPI_MM_AMN
   integer(c_int), bind(c, name='mf_mpi_ww_sum_f'),     protected :: MPI_WW_SUM
   integer(c_int), bind(c, name='mf_mpi_ww_amx_f'),     protected :: MPI_WW_AMX
   integer(c_int), bind(c, name='mf_mpi_ww_amn_f'),     protected :: MPI_WW_AMN

   integer(c_int), bind(c, name='mf_mpi_dd_sum_f'),     protected :: MPI_DD_SUM
   integer(c_int), bind(c, name='mf_mpi_dd_amx_f'),     protected :: MPI_DD_AMX
   integer(c_int), bind(c, name='mf_mpi_dd_amn_f'),     protected :: MPI_DD_AMN
   integer(c_int), bind(c, name='mf_mpi_zz_sum_f'),     protected :: MPI_ZZ_SUM
   integer(c_int), bind(c, name='mf_mpi_zz_amx_f'),     protected :: MPI_ZZ_AMX
   integer(c_int), bind(c, name='mf_mpi_zz_amn_f'),     protected :: MPI_ZZ_AMN

   interface
      subroutine multifloats_mpi_init() bind(c, name='multifloats_mpi_init')
      end subroutine multifloats_mpi_init
   end interface
end module multifloats_mpi_f
