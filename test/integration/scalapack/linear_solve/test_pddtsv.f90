program test_pddtsv
    ! ScaLAPACK 1D-distributed tridiagonal solver. Diagonals (DL, D,
    ! DU) are distributed in column-block fashion across a 1xNPROCS
    ! BLACS grid; B is held in P-by-1 layout. Reference is the serial
    ! quad LAPACK dgtsv.
    !
    ! Each rank generates its local DL/D/DU slab from the global vectors
    ! at REAL(KIND=16). Diagonal D is boosted by 2*N to keep the system
    ! diagonally dominant (a precondition for pddtsv with no pivoting).
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgtsv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_nproc, &
                                my1d_context, my1d_npcol, my1d_col, &
                                descinit_1d
    use test_data,        only: gen_vector_quad, gen_matrix_quad
    use target_scalapack, only: target_name, target_eps, target_pddtsv
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer :: i, n, nb, info, info_ref, lwork, ig, jl, lldB, locm_b
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer, parameter :: tag = 7611
    integer :: desca(9), descb(9)
    real(ep), allocatable :: dl_glob(:), d_glob(:), du_glob(:), B_glob(:,:)
    real(ep), allocatable :: dl_loc(:),  d_loc(:),  du_loc(:),  B_loc(:,:)
    real(ep), allocatable :: B_got(:,:), B_ref(:,:), buf(:,:)
    real(ep), allocatable :: dl_ref(:), d_ref(:), du_ref(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pddtsv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        ! Choose a block size that distributes evenly: NB = N / NPROCS
        ! (rounded up if needed). Restriction: P*NB >= N.
        nb = (n + my_nproc - 1) / my_nproc

        ! Generate global diagonals. Boost D for diagonal dominance.
        ! Convention split: pddtsv uses DL(2..N) for the sub-diagonal
        ! (DL(1) unused); dgtsv uses DL(1..N-1) for the same physical
        ! entries. Generate a single physical sub-diagonal and lay it
        ! out for pddtsv; the reference path shifts it back at call.
        call gen_vector_quad(n,     dl_glob, seed = 18001 + 31*i)
        call gen_vector_quad(n,     d_glob,  seed = 18011 + 31*i)
        call gen_vector_quad(n,     du_glob, seed = 18021 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 18031 + 31*i)
        d_glob = d_glob + real(2 * n, ep)
        ! Zero the unused slots so the two conventions agree on
        ! "off-bandwidth" elements.
        dl_glob(1) = 0.0_ep
        du_glob(n) = 0.0_ep

        ! Local slabs of size NB per rank — last rank may have fewer
        ! valid entries but the array length is always NB.
        allocate(dl_loc(nb), d_loc(nb), du_loc(nb))
        dl_loc = 0.0_ep; d_loc = 0.0_ep; du_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb
                owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    dl_loc(jl) = dl_glob(ig)
                    d_loc(jl)  = d_glob(ig)
                    du_loc(jl) = du_glob(ig)
                end if
            end block
        end do

        ! B is held in NB-row blocks per process.
        lldB = nb
        locm_b = nb
        allocate(B_loc(nb, nrhs))
        B_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb
                owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    B_loc(jl, 1:nrhs) = B_glob(ig, 1:nrhs)
                end if
            end block
        end do

        ! A uses DTYPE=501 (1xP); B must use DTYPE=502 (Px1) per the
        ! pddtsv documentation.
        call descinit_1d(desca, 501, n, nb, 0, my1d_context, nb,   info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        ! pddtsv documented LWORK >= 12*NPCOL + 3*NB + max(10*NPCOL+4*NRHS, 8*NPCOL).
        lwork = 12 * my1d_npcol + 3 * nb + max(10 * my1d_npcol + 4 * nrhs, 8 * my1d_npcol)
        allocate(work(max(1, lwork)))

        call target_pddtsv(n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                           B_loc, 1, descb, work, lwork, info)

        ! Gather B_loc → B_got on rank 0.
        bytes_per_elem = storage_size(0.0_ep) / 8
        if (my_rank == 0) then
            allocate(B_got(n, nrhs))
            B_got = 0.0_ep
            do src_rank = 0, my_nproc - 1
                allocate(buf(nb, nrhs))
                if (src_rank == my_rank) then
                    buf = B_loc
                else
                    call mpi_recv(buf, nb * nrhs * bytes_per_elem, mpi_byte, &
                                  src_rank, tag, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    block
                        integer :: blk, owner
                        blk = (ig - 1) / nb
                        owner = blk
                        if (owner == src_rank) then
                            jl = mod(ig - 1, nb) + 1
                            B_got(ig, 1:nrhs) = buf(jl, 1:nrhs)
                        end if
                    end block
                end do
                deallocate(buf)
            end do
        else
            peer = 0
            call mpi_send(B_loc, nb * nrhs * bytes_per_elem, mpi_byte, &
                          peer, tag, mpi_comm_world, ierr)
        end if

        if (my_rank == 0) then
            allocate(dl_ref(n), d_ref(n), du_ref(n), B_ref(n, nrhs))
            ! dgtsv expects sub-diagonal in dl_ref(1..n-1); pddtsv put
            ! it in dl_glob(2..n). Shift on copy. du is the same in
            ! both conventions.
            dl_ref(1:n-1) = dl_glob(2:n); dl_ref(n) = 0.0_ep
            d_ref = d_glob; du_ref = du_glob
            B_ref = B_glob
            call dgtsv(n, nrhs, dl_ref, d_ref, du_ref, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(dl_ref, d_ref, du_ref, B_ref, B_got)
        end if
        deallocate(dl_loc, d_loc, du_loc, B_loc, work, &
                   dl_glob, d_glob, du_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pddtsv
