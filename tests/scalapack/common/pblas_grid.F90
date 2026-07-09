! BLACS grid setup and block-cyclic index helpers for PBLAS tests.
!
! Initializes MPI + BLACS once per test program (grid_init) and tears
! down at exit (grid_exit). Process grid is 2D, shaped as close to
! square as possible — a 2x2 for mpirun -np 4, 1x1 for -np 1, 1x3 for
! prime np, etc. Row-major ordering matches BLACS_GRIDINIT('R').
!
! numroc_local and descinit_local are pure Fortran implementations.
! They duplicate a few lines from ScaLAPACK's TOOLS/numroc.f and
! TOOLS/descinit.f, but inlining them here avoids relying on the
! migrated ScaLAPACK library being present (tests/pblas only depends
! on ${LIB_PREFIX}blacs + ${LIB_PREFIX}pblas, not scalapack).

module pblas_grid
    use mpi
    implicit none
    private

    public :: grid_init, grid_exit
    public :: my_context, my_nprow, my_npcol, my_row, my_col
    public :: my_rank, my_nproc
    public :: pick_grid_shape
    public :: numroc_local
    public :: descinit_local
    public :: g2l
    ! Companion 1xNPROCS BLACS context used by the ScaLAPACK 1D-
    ! distributed routines (the tridiagonal / banded solver families).
    ! The CTXT-handle and grid info are exposed as my1d_context / my1d_*
    ! and the descriptor is built via descinit_1d. Ranks live on row 0
    ! of the 1D grid so my1d_col == my_rank.
    public :: my1d_context, my1d_npcol, my1d_col
    public :: descinit_1d

    integer, save :: my_context = -1
    integer, save :: my_nprow   = 0
    integer, save :: my_npcol   = 0
    integer, save :: my_row     = -1
    integer, save :: my_col     = -1
    integer, save :: my_rank    = -1
    integer, save :: my_nproc   = 0
    integer, save :: my1d_context = -1
    integer, save :: my1d_npcol   = 0
    integer, save :: my1d_col     = -1

    interface
        subroutine blacs_pinfo(mypnum, nprocs)
            integer, intent(out) :: mypnum, nprocs
        end subroutine
        subroutine blacs_get(icontxt, what, val)
            integer, intent(in)  :: icontxt, what
            integer, intent(out) :: val
        end subroutine
        subroutine blacs_gridinit(icontxt, order, nprow, npcol)
            integer,      intent(inout) :: icontxt
            character(*), intent(in)    :: order
            integer,      intent(in)    :: nprow, npcol
        end subroutine
        subroutine blacs_gridinfo(icontxt, nprow, npcol, myrow, mycol)
            integer, intent(in)  :: icontxt
            integer, intent(out) :: nprow, npcol, myrow, mycol
        end subroutine
        subroutine blacs_gridexit(icontxt)
            integer, intent(in) :: icontxt
        end subroutine
        subroutine blacs_exit(continue_)
            integer, intent(in) :: continue_
        end subroutine
    end interface

contains

    subroutine grid_init()
        integer :: ierr, ctxt0, ctxt1d, dummy_nprow, dummy_row
        logical :: mpi_started

        call mpi_initialized(mpi_started, ierr)
        if (.not. mpi_started) call mpi_init(ierr)

        call blacs_pinfo(my_rank, my_nproc)
        if (my_nproc < 1) call mpi_comm_size(mpi_comm_world, my_nproc, ierr)

        call pick_grid_shape(my_nproc, my_nprow, my_npcol)

        ! BLACS_GET with what=0 retrieves the default system context
        ! (derived from MPI_COMM_WORLD). Reusing it as the handle we
        ! init the grid on keeps everything under one communicator.
        call blacs_get(-1, 0, ctxt0)
        my_context = ctxt0
        call blacs_gridinit(my_context, 'R', my_nprow, my_npcol)
        call blacs_gridinfo(my_context, my_nprow, my_npcol, my_row, my_col)

        ! Companion 1xNPROCS context for the ScaLAPACK 1D-distribution
        ! routines (tridiagonal / banded solver families). Same MPI
        ! ranks, just laid out as a single row of NPROCS columns.
        call blacs_get(-1, 0, ctxt1d)
        my1d_context = ctxt1d
        call blacs_gridinit(my1d_context, 'R', 1, my_nproc)
        call blacs_gridinfo(my1d_context, dummy_nprow, my1d_npcol, &
                            dummy_row, my1d_col)
    end subroutine grid_init

    subroutine grid_exit()
        if (my1d_col >= 0) then
            call blacs_gridexit(my1d_context)
        end if
        if (my_row >= 0 .and. my_col >= 0) then
            call blacs_gridexit(my_context)
        end if
        ! continue_=1 leaves MPI alive for MPI_Finalize to tear down;
        ! continue_=0 would finalize MPI inside BLACS.
        call blacs_exit(1)
        call mpi_finalize_if_needed()
    end subroutine grid_exit

    subroutine mpi_finalize_if_needed()
        integer :: ierr
        logical :: finalized
        call mpi_finalized(finalized, ierr)
        if (.not. finalized) call mpi_finalize(ierr)
    end subroutine

    ! Factor `np` into (nprow, npcol) with nprow as large as possible
    ! while staying <= sqrt(np). Produces a square grid for perfect
    ! squares, a near-square for non-squares, and 1×np for primes.
    subroutine pick_grid_shape(np, nprow, npcol)
        integer, intent(in)  :: np
        integer, intent(out) :: nprow, npcol
        integer :: i
        nprow = 1
        do i = 2, int(sqrt(real(np))) + 1
            if (mod(np, i) == 0 .and. i * i <= np) nprow = i
        end do
        npcol = np / nprow
    end subroutine pick_grid_shape

    ! Block-cyclic row/column count for a given process.
    !   n       : global dimension
    !   nb      : block size
    !   iproc   : this process's row or column index in the grid
    !   isrcproc: process holding the first block (normally 0)
    !   nprocs  : number of processes in that dimension
    pure function numroc_local(n, nb, iproc, isrcproc, nprocs) result(locn)
        integer, intent(in) :: n, nb, iproc, isrcproc, nprocs
        integer :: locn, mydist, nblocks, extrablks
        mydist   = mod(nprocs + iproc - isrcproc, nprocs)
        nblocks  = n / nb
        locn     = (nblocks / nprocs) * nb
        extrablks = mod(nblocks, nprocs)
        if (mydist < extrablks) then
            locn = locn + nb
        else if (mydist == extrablks) then
            locn = locn + mod(n, nb)
        end if
    end function numroc_local

    ! Build a ScaLAPACK descriptor (9 integers) for a distributed
    ! matrix. Matches TOOLS/descinit.f byte-for-byte for the fields
    ! we care about; the only validation is lld >= locm on the owning
    ! process — PBLAS will catch anything else.
    subroutine descinit_local(desc, m, n, mb, nb, irsrc, icsrc, &
                              ctxt, lld, info)
        integer, intent(out) :: desc(9)
        integer, intent(in)  :: m, n, mb, nb, irsrc, icsrc, ctxt, lld
        integer, intent(out) :: info
        info     = 0
        desc(1)  = 1        ! DTYPE_A = 1 (dense matrix)
        desc(2)  = ctxt     ! context
        desc(3)  = m        ! global rows
        desc(4)  = n        ! global cols
        desc(5)  = mb       ! row block size
        desc(6)  = nb       ! column block size
        desc(7)  = irsrc    ! process row that owns A(1,1)
        desc(8)  = icsrc    ! process col that owns A(1,1)
        desc(9)  = max(1, lld)
    end subroutine descinit_local

    ! Build a 1D ScaLAPACK descriptor (7 integers — sized 9 to share
    ! the same array shape as 2D contexts).
    !   dtype : 501 (1xP) for A in tridiagonal/banded routines;
    !           502 (Px1) is the only valid choice for B.
    !   n     : global problem size
    !   nb    : block size
    !   isrc  : source process column (normally 0)
    !   ctxt  : BLACS context (must refer to a 1D grid)
    !   lld   : local leading dimension (used for B; ignored for A)
    subroutine descinit_1d(desc, dtype, n, nb, isrc, ctxt, lld, info)
        integer, intent(out) :: desc(9)
        integer, intent(in)  :: dtype, n, nb, isrc, ctxt, lld
        integer, intent(out) :: info
        info    = 0
        desc(1) = dtype     ! DTYPE_A: 501 (1xP) or 502 (Px1)
        desc(2) = ctxt      ! 1D context
        desc(3) = n         ! global size
        desc(4) = nb        ! block size
        desc(5) = isrc      ! source process index
        desc(6) = max(1, lld)
        desc(7) = 0
        desc(8) = 0
        desc(9) = 0
    end subroutine descinit_1d

    ! Block-cyclic global-to-local mapping. Given a global 1-based
    ! index, block size, and number of processes in that dimension,
    ! returns the owning process index and the local 1-based index.
    pure subroutine g2l(ig, b, p, owner, il)
        integer, intent(in)  :: ig, b, p
        integer, intent(out) :: owner, il
        integer :: zero_based, block_idx
        zero_based = ig - 1
        block_idx  = zero_based / b
        owner      = mod(block_idx, p)
        il         = (block_idx / p) * b + mod(zero_based, b) + 1
    end subroutine g2l

end module pblas_grid
