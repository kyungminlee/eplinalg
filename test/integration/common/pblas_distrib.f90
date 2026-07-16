! Block-cyclic distribution helpers for PBLAS tests.
!
! The tests follow a mirrored-global pattern: every rank independently
! materializes the full global vector/matrix at REAL(KIND=ep) using
! the same seed (so they have identical data), then copies its local
! block-cyclic piece into a smaller local array that PBLAS operates on.
! After the PBLAS call, local results are gathered back to rank 0 for
! comparison against serial refblas_quad applied to the original
! global.
!
! Raw MPI_BYTE transfers are used for the gather because MPI_REAL16 is
! not a portable predefined type across MPI implementations. A local
! element of REAL(KIND=16) is 16 bytes, COMPLEX(KIND=16) is 32 bytes.
!
! The scatter/gather layout assumes block-cyclic over a 2D grid with
! process (0,0) owning A(1,1), block sizes mb × nb, and column-major
! local storage (leading dimension = local row count).

module pblas_distrib
    use prec_kinds, only: ep
    use test_data,  only: gen_vector_quad, gen_matrix_quad, &
                          gen_vector_complex, gen_matrix_complex
    use pblas_grid, only: my_nprow, my_npcol, my_row, my_col, &
                          numroc_local, g2l
    use mpi
    implicit none
    private

    public :: gen_distrib_vector, gen_distrib_matrix, gen_distrib_vector_row
    public :: gather_vector, gather_matrix, gather_vector_row
    public :: gen_distrib_vector_z, gen_distrib_matrix_z
    public :: gather_vector_z, gather_matrix_z

contains

    ! ── Real vector ─────────────────────────────────────────────────
    ! "Column-vector" layout: distributed as a single column of an
    ! n × 1 matrix. Rows are split across `nprow` processes in `mb`
    ! sized chunks; only column 0 holds data.
    subroutine gen_distrib_vector(n, mb, x_loc, x_glob, seed)
        integer,  intent(in) :: n, mb, seed
        real(ep), intent(out), allocatable :: x_loc(:), x_glob(:)
        integer :: loc_n, ig, owner_r, il

        call gen_vector_quad(n, x_glob, seed)
        ! Allocate buffer of size LLD on every rank — the descriptor's
        ! LLD = max(1, numroc(n, mb, my_row, 0, my_nprow)) doesn't
        ! depend on my_col, so wrappers like target_pdcopy reading
        ! x(1:lld*nb) on non-column-0 ranks would otherwise overrun
        ! a smaller `max(1, 0)` placeholder. Detected by valgrind under
        ! real distributed MPI; singleton MPI never tripped it because
        ! my_nprow=my_npcol=1 made the my_col gate inert.
        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        allocate(x_loc(max(1, loc_n)))
        x_loc = 0.0_ep
        if (my_col == 0 .and. loc_n > 0) then
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                if (owner_r == my_row) x_loc(il) = x_glob(ig)
            end do
        end if
    end subroutine gen_distrib_vector

    ! "Row-vector" layout: 1 × n matrix. Columns distributed across
    ! `npcol` processes in `nb`-sized chunks; only row 0 holds data.
    subroutine gen_distrib_vector_row(n, nb, x_loc, x_glob, seed)
        integer,  intent(in) :: n, nb, seed
        real(ep), intent(out), allocatable :: x_loc(:), x_glob(:)
        integer :: loc_n, ig, owner_c, jl

        call gen_vector_quad(n, x_glob, seed)
        ! Allocate full LLD-size buffer regardless of my_row — see the
        ! note in gen_distrib_vector for the rationale (descriptor's
        ! LLD doesn't depend on my_row, so wrappers reading the buffer
        ! on non-row-0 ranks must see the correct size).
        loc_n = numroc_local(n, nb, my_col, 0, my_npcol)
        allocate(x_loc(max(1, loc_n)))
        x_loc = 0.0_ep
        if (my_row == 0 .and. loc_n > 0) then
            do ig = 1, n
                call g2l(ig, nb, my_npcol, owner_c, jl)
                if (owner_c == my_col) x_loc(jl) = x_glob(ig)
            end do
        end if
    end subroutine gen_distrib_vector_row

    ! ── Real matrix ─────────────────────────────────────────────────
    subroutine gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed)
        integer,  intent(in) :: m, n, mb, nb, seed
        real(ep), intent(out), allocatable :: A_loc(:,:), A_glob(:,:)
        integer :: locm, locn, ig, jg, owner_r, owner_c, il, jl

        call gen_matrix_quad(m, n, A_glob, seed)
        locm = numroc_local(m, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol)
        allocate(A_loc(max(1, locm), max(1, locn)))
        A_loc = 0.0_ep
        if (locm > 0 .and. locn > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, m
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_glob(ig, jg)
                end do
            end do
        end if
    end subroutine gen_distrib_matrix

    ! ── Complex vector / matrix ─────────────────────────────────────
    subroutine gen_distrib_vector_z(n, mb, x_loc, x_glob, seed)
        integer,     intent(in) :: n, mb, seed
        complex(ep), intent(out), allocatable :: x_loc(:), x_glob(:)
        integer :: loc_n, ig, owner_r, il

        call gen_vector_complex(n, x_glob, seed)
        ! Allocate full LLD-size buffer regardless of my_col — see
        ! gen_distrib_vector for the rationale.
        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        allocate(x_loc(max(1, loc_n)))
        x_loc = (0.0_ep, 0.0_ep)
        if (my_col == 0 .and. loc_n > 0) then
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                if (owner_r == my_row) x_loc(il) = x_glob(ig)
            end do
        end if
    end subroutine gen_distrib_vector_z

    subroutine gen_distrib_matrix_z(m, n, mb, nb, A_loc, A_glob, seed)
        integer,     intent(in) :: m, n, mb, nb, seed
        complex(ep), intent(out), allocatable :: A_loc(:,:), A_glob(:,:)
        integer :: locm, locn, ig, jg, owner_r, owner_c, il, jl

        call gen_matrix_complex(m, n, A_glob, seed)
        locm = numroc_local(m, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol)
        allocate(A_loc(max(1, locm), max(1, locn)))
        A_loc = (0.0_ep, 0.0_ep)
        if (locm > 0 .and. locn > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, m
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_glob(ig, jg)
                end do
            end do
        end if
    end subroutine gen_distrib_matrix_z

    ! ── Gather: local block-cyclic → rank-0 global ──────────────────
    ! Sends each (prow, pcol)'s local slab to rank 0 as raw bytes;
    ! rank 0 deposits it into the correct global positions using the
    ! block-cyclic map. Non-zero ranks allocate nothing on x_glob.
    subroutine gather_vector(n, mb, x_loc, x_glob)
        integer,  intent(in) :: n, mb
        real(ep), intent(in) :: x_loc(:)
        real(ep), intent(out), allocatable :: x_glob(:)
        integer :: ig, il, owner_r, ierr, peer, bytes_per_elem
        integer :: src_rank, loc_n, pr
        integer, parameter :: tag = 6411
        real(ep), allocatable :: buf(:)

        bytes_per_elem = storage_size(x_loc(1)) / 8

        if (my_rank_global() == 0) then
            allocate(x_glob(n))
            x_glob = 0.0_ep
            do pr = 0, my_nprow - 1
                loc_n = numroc_local(n, mb, pr, 0, my_nprow)
                if (loc_n == 0) cycle
                allocate(buf(loc_n))
                if (pr == my_row .and. my_col == 0) then
                    buf = x_loc(1:loc_n)
                else
                    src_rank = pr * my_npcol  ! (pr, 0) in row-major
                    call mpi_recv(buf, loc_n * bytes_per_elem, mpi_byte, &
                                  src_rank, tag, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == pr) x_glob(ig) = buf(il)
                end do
                deallocate(buf)
            end do
        else
            if (my_col == 0) then
                loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
                if (loc_n > 0) then
                    peer = 0
                    call mpi_send(x_loc, loc_n * bytes_per_elem, mpi_byte, &
                                  peer, tag, mpi_comm_world, ierr)
                end if
            end if
        end if
    end subroutine gather_vector

    subroutine gather_vector_row(n, nb, x_loc, x_glob)
        integer,  intent(in) :: n, nb
        real(ep), intent(in) :: x_loc(:)
        real(ep), intent(out), allocatable :: x_glob(:)
        integer :: ig, jl, owner_c, ierr, peer, bytes_per_elem
        integer :: src_rank, loc_n, pc
        integer, parameter :: tag = 6412
        real(ep), allocatable :: buf(:)

        bytes_per_elem = storage_size(x_loc(1)) / 8

        if (my_rank_global() == 0) then
            allocate(x_glob(n))
            x_glob = 0.0_ep
            do pc = 0, my_npcol - 1
                loc_n = numroc_local(n, nb, pc, 0, my_npcol)
                if (loc_n == 0) cycle
                allocate(buf(loc_n))
                if (pc == my_col .and. my_row == 0) then
                    buf = x_loc(1:loc_n)
                else
                    src_rank = pc  ! (0, pc) in row-major
                    call mpi_recv(buf, loc_n * bytes_per_elem, mpi_byte, &
                                  src_rank, tag, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    call g2l(ig, nb, my_npcol, owner_c, jl)
                    if (owner_c == pc) x_glob(ig) = buf(jl)
                end do
                deallocate(buf)
            end do
        else
            if (my_row == 0) then
                loc_n = numroc_local(n, nb, my_col, 0, my_npcol)
                if (loc_n > 0) then
                    peer = 0
                    call mpi_send(x_loc, loc_n * bytes_per_elem, mpi_byte, &
                                  peer, tag, mpi_comm_world, ierr)
                end if
            end if
        end if
    end subroutine gather_vector_row

    subroutine gather_matrix(m, n, mb, nb, A_loc, A_glob)
        integer,  intent(in) :: m, n, mb, nb
        real(ep), intent(in) :: A_loc(:,:)
        real(ep), intent(out), allocatable :: A_glob(:,:)
        integer :: ig, jg, il, jl, owner_r, owner_c, ierr, peer
        integer :: bytes_per_elem, src_rank, locm, locn, pr, pc
        integer, parameter :: tag = 6413
        real(ep), allocatable :: buf(:,:)

        bytes_per_elem = storage_size(A_loc(1,1)) / 8

        if (my_rank_global() == 0) then
            allocate(A_glob(m, n))
            A_glob = 0.0_ep
            do pr = 0, my_nprow - 1
                do pc = 0, my_npcol - 1
                    locm = numroc_local(m, mb, pr, 0, my_nprow)
                    locn = numroc_local(n, nb, pc, 0, my_npcol)
                    if (locm == 0 .or. locn == 0) cycle
                    allocate(buf(locm, locn))
                    if (pr == my_row .and. pc == my_col) then
                        buf = A_loc(1:locm, 1:locn)
                    else
                        src_rank = pr * my_npcol + pc
                        call mpi_recv(buf, locm * locn * bytes_per_elem, &
                                      mpi_byte, src_rank, tag, &
                                      mpi_comm_world, mpi_status_ignore, ierr)
                    end if
                    do jg = 1, n
                        call g2l(jg, nb, my_npcol, owner_c, jl)
                        if (owner_c /= pc) cycle
                        do ig = 1, m
                            call g2l(ig, mb, my_nprow, owner_r, il)
                            if (owner_r == pr) A_glob(ig, jg) = buf(il, jl)
                        end do
                    end do
                    deallocate(buf)
                end do
            end do
        else
            locm = numroc_local(m, mb, my_row, 0, my_nprow)
            locn = numroc_local(n, nb, my_col, 0, my_npcol)
            if (locm > 0 .and. locn > 0) then
                peer = 0
                call mpi_send(A_loc, locm * locn * bytes_per_elem, mpi_byte, &
                              peer, tag, mpi_comm_world, ierr)
            end if
        end if
    end subroutine gather_matrix

    ! Complex gather variants — same pattern, 2× bytes per element.
    subroutine gather_vector_z(n, mb, x_loc, x_glob)
        integer,     intent(in) :: n, mb
        complex(ep), intent(in) :: x_loc(:)
        complex(ep), intent(out), allocatable :: x_glob(:)
        integer :: ig, il, owner_r, ierr, peer, bytes_per_elem
        integer :: src_rank, loc_n, pr
        integer, parameter :: tag = 6421
        complex(ep), allocatable :: buf(:)

        bytes_per_elem = storage_size(x_loc(1)) / 8

        if (my_rank_global() == 0) then
            allocate(x_glob(n))
            x_glob = (0.0_ep, 0.0_ep)
            do pr = 0, my_nprow - 1
                loc_n = numroc_local(n, mb, pr, 0, my_nprow)
                if (loc_n == 0) cycle
                allocate(buf(loc_n))
                if (pr == my_row .and. my_col == 0) then
                    buf = x_loc(1:loc_n)
                else
                    src_rank = pr * my_npcol
                    call mpi_recv(buf, loc_n * bytes_per_elem, mpi_byte, &
                                  src_rank, tag, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == pr) x_glob(ig) = buf(il)
                end do
                deallocate(buf)
            end do
        else
            if (my_col == 0) then
                loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
                if (loc_n > 0) then
                    peer = 0
                    call mpi_send(x_loc, loc_n * bytes_per_elem, mpi_byte, &
                                  peer, tag, mpi_comm_world, ierr)
                end if
            end if
        end if
    end subroutine gather_vector_z

    subroutine gather_matrix_z(m, n, mb, nb, A_loc, A_glob)
        integer,     intent(in) :: m, n, mb, nb
        complex(ep), intent(in) :: A_loc(:,:)
        complex(ep), intent(out), allocatable :: A_glob(:,:)
        integer :: ig, jg, il, jl, owner_r, owner_c, ierr, peer
        integer :: bytes_per_elem, src_rank, locm, locn, pr, pc
        integer, parameter :: tag = 6423
        complex(ep), allocatable :: buf(:,:)

        bytes_per_elem = storage_size(A_loc(1,1)) / 8

        if (my_rank_global() == 0) then
            allocate(A_glob(m, n))
            A_glob = (0.0_ep, 0.0_ep)
            do pr = 0, my_nprow - 1
                do pc = 0, my_npcol - 1
                    locm = numroc_local(m, mb, pr, 0, my_nprow)
                    locn = numroc_local(n, nb, pc, 0, my_npcol)
                    if (locm == 0 .or. locn == 0) cycle
                    allocate(buf(locm, locn))
                    if (pr == my_row .and. pc == my_col) then
                        buf = A_loc(1:locm, 1:locn)
                    else
                        src_rank = pr * my_npcol + pc
                        call mpi_recv(buf, locm * locn * bytes_per_elem, &
                                      mpi_byte, src_rank, tag, &
                                      mpi_comm_world, mpi_status_ignore, ierr)
                    end if
                    do jg = 1, n
                        call g2l(jg, nb, my_npcol, owner_c, jl)
                        if (owner_c /= pc) cycle
                        do ig = 1, m
                            call g2l(ig, mb, my_nprow, owner_r, il)
                            if (owner_r == pr) A_glob(ig, jg) = buf(il, jl)
                        end do
                    end do
                    deallocate(buf)
                end do
            end do
        else
            locm = numroc_local(m, mb, my_row, 0, my_nprow)
            locn = numroc_local(n, nb, my_col, 0, my_npcol)
            if (locm > 0 .and. locn > 0) then
                peer = 0
                call mpi_send(A_loc, locm * locn * bytes_per_elem, mpi_byte, &
                              peer, tag, mpi_comm_world, ierr)
            end if
        end if
    end subroutine gather_matrix_z

    ! Helper: MPI_COMM_WORLD rank. BLACS grid's own rank mapping is
    ! row-major (rank = my_row * my_npcol + my_col), which also equals
    ! the MPI rank when BLACS gets its context from MPI_COMM_WORLD.
    integer function my_rank_global()
        integer :: r, ierr
        call mpi_comm_rank(mpi_comm_world, r, ierr)
        my_rank_global = r
    end function my_rank_global

end module pblas_distrib
