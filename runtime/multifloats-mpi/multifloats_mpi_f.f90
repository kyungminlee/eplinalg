! multifloats_mpi_f -- Fortran-side MPI handles for the real-precision
! multifloats datatype and the real + complex reduction operators
! registered by multifloats_mpi_init().
!
! Migrated multifloats Fortran code (e.g. mmumps) calls MPI directly
! with statements like ``CALL MPI_SEND(..., MPI_FLOAT64X2, ...)`` and
! ``CALL MPI_REDUCE(..., MPI_DD_SUM, ...)``. The handles are created in
! C++ at runtime, not by mpif.h, so Fortran needs an explicit module
! that exposes them. Each public name below is a default-INTEGER
! (matching MPI_Fint on Linux/x86_64) bound to the C symbol
! multifloats_mpi.cpp populates via MPI_Type_c2f / MPI_Op_c2f. Values
! are 0 until ``multifloats_mpi_init`` runs — that happens
! automatically the first time BLACS bootstraps MPI, and can be
! triggered explicitly via ``CALL multifloats_mpi_init`` for consumers
! (like MUMPS) that don't go through BLACS.
!
module multifloats_mpi_f
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   private

   public :: MPI_FLOAT64X2, MPI_COMPLEX64X2
   public :: MPI_DD_SUM, MPI_DD_AMX, MPI_DD_AMN
   public :: MPI_ZZ_SUM, MPI_ZZ_AMX, MPI_ZZ_AMN
   public :: multifloats_mpi_init

   integer(c_int), bind(c, name='mf_mpi_float64x2_f'),  protected :: MPI_FLOAT64X2
   integer(c_int), bind(c, name='mf_mpi_complex64x2_f'),protected :: MPI_COMPLEX64X2
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
