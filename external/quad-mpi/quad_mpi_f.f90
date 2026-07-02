! quad_mpi_f -- Fortran-side handles for the custom MPI reduction ops
! used by the migrated kind16 (REAL(KIND=16) / __float128) MUMPS.
!
! The datatypes themselves stay standard (MPI_REAL16 / MPI_COMPLEX32,
! from mpif.h) -- only the reduction *ops* are custom, because Intel
! MPI's builtin MPI_SUM / MPI_MAX / MPI_MIN have no 16-byte-real kernel.
! Migrated qmumps / xmumps code calls MPI directly with statements like
! ``CALL MPI_REDUCE(..., MPI_REAL16, MPI_QQ_SUM, ...)``; the ops are
! created in C at runtime (not by mpif.h), so Fortran needs this module
! to expose them. Each public name is a default-INTEGER (matching
! MPI_Fint on Linux/x86_64) bound to the C symbol quad_mpi.c populates
! via MPI_Op_c2f. Values are 0 until ``quad_mpi_init`` runs -- that
! happens automatically the first time BLACS bootstraps MPI, and can be
! triggered explicitly via ``CALL quad_mpi_init`` for consumers (like
! MUMPS) that reach a reduction before any BLACS call.
!
module quad_mpi_f
   use, intrinsic :: iso_c_binding, only: c_int
   implicit none
   private

   public :: MPI_QQ_SUM, MPI_QQ_AMX, MPI_QQ_AMN
   public :: MPI_XX_SUM, MPI_XX_AMX, MPI_XX_AMN
   public :: quad_mpi_init

   integer(c_int), bind(c, name='q_mpi_qq_sum_f'), protected :: MPI_QQ_SUM
   integer(c_int), bind(c, name='q_mpi_qq_amx_f'), protected :: MPI_QQ_AMX
   integer(c_int), bind(c, name='q_mpi_qq_amn_f'), protected :: MPI_QQ_AMN
   integer(c_int), bind(c, name='q_mpi_xx_sum_f'), protected :: MPI_XX_SUM
   integer(c_int), bind(c, name='q_mpi_xx_amx_f'), protected :: MPI_XX_AMX
   integer(c_int), bind(c, name='q_mpi_xx_amn_f'), protected :: MPI_XX_AMN

   interface
      subroutine quad_mpi_init() bind(c, name='quad_mpi_init')
      end subroutine quad_mpi_init
   end interface
end module quad_mpi_f
