program test_pdpttrsv
    ! Validate the SPD-tridiagonal half-solves by composing them to recover
    ! the full pdpttrs result. PDPTTRF stores A = L*D*L^T with the LDL^T
    ! diagonal D' in the modified d_loc and L's bidiagonal subdiagonal
    ! folded into AF; PDPTTRSV(UPLO='L') applies L^-1, PDPTTRSV(UPLO='U')
    ! applies L^-T. The full pdpttrs solve corresponds to
    !     X = L^-T * (D'^-1 .* (L^-1 * B))
    ! so the per-rank entrywise divide by d_loc(jl) belongs between the
    ! two halves. (The naïve trsv(L) ∘ trsv(U) without that divide
    ! disagrees by 60-200x — that's what the prior smoke test got wrong.)
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_npcol, my1d_col, descinit_1d
    use test_data,         only: gen_vector_quad, gen_matrix_quad
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pdpttrf, target_pdpttrs, target_pdpttrsv
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: tag = 8801
    integer :: i, n, nb, info, lwork, laf, ig, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer :: desca(9), descb(9)
    real(ep), allocatable :: d_glob(:), e_glob(:), B_glob(:,:)
    real(ep), allocatable :: d_loc(:),  e_loc(:)
    real(ep), allocatable :: B_full(:,:), B_half(:,:)
    real(ep), allocatable :: af(:), work(:)
    real(ep), allocatable :: B_full_got(:,:), B_half_got(:,:), buf(:,:)
    real(ep) :: err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pdpttrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, d_glob, seed = 33301 + 31*i)
        call gen_vector_quad(n, e_glob, seed = 33311 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 33321 + 31*i)
        d_glob = abs(d_glob) + real(2 * n, ep)
        e_glob = e_glob * 0.05_ep
        e_glob(n) = 0.0_ep

        allocate(d_loc(nb), e_loc(nb))
        d_loc = 0.0_ep; e_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    d_loc(jl) = d_glob(ig)
                    e_loc(jl) = e_glob(ig)
                end if
            end block
        end do

        lldB = nb
        allocate(B_full(nb, nrhs), B_half(nb, nrhs))
        B_full = 0.0_ep; B_half = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    B_full(jl, 1:nrhs) = B_glob(ig, 1:nrhs)
                    B_half(jl, 1:nrhs) = B_glob(ig, 1:nrhs)
                end if
            end block
        end do

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, nb,   info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        laf = 12 * my1d_npcol + 3 * nb
        allocate(af(laf))

        call target_pdpttrf(n, d_loc, e_loc, 1, desca, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrf(n, d_loc, e_loc, 1, desca, af, laf, work, lwork, info)
        deallocate(work)

        ! Full solve (oracle).
        call target_pdpttrs(n, nrhs, d_loc, e_loc, 1, desca, &
                            B_full, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrs(n, nrhs, d_loc, e_loc, 1, desca, &
                            B_full, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        ! Half-solve composition: trsv(L) -> divide by D' -> trsv(U).
        call target_pdpttrsv('L', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_half, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpttrsv('L', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        ! Mid-step divide between trsv halves. Mirrors upstream
        ! pdpttrs.f:705-716: ranks 0..NPCOL-2 divide entries 1..ODD_SIZE
        ! (= nb-1) by d_loc, then divide the boundary entry nb by the
        ! reduced-system pivot stored at AF(nb+1). The last rank
        ! divides all nb entries by d_loc (no boundary).
        if (my1d_col < my1d_npcol - 1) then
            do jl = 1, nb - 1
                if (d_loc(jl) /= 0.0_ep) &
                    B_half(jl, 1:nrhs) = B_half(jl, 1:nrhs) / d_loc(jl)
            end do
            if (af(nb + 1) /= 0.0_ep) &
                B_half(nb, 1:nrhs) = B_half(nb, 1:nrhs) / af(nb + 1)
        else
            do jl = 1, nb
                if (d_loc(jl) /= 0.0_ep) &
                    B_half(jl, 1:nrhs) = B_half(jl, 1:nrhs) / d_loc(jl)
            end do
        end if

        allocate(work(1))
        call target_pdpttrsv('U', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdpttrsv('U', n, nrhs, d_loc, e_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)

        bytes_per_elem = storage_size(0.0_ep) / 8
        if (my_rank == 0) then
            allocate(B_full_got(n, nrhs), B_half_got(n, nrhs))
            B_full_got = 0.0_ep; B_half_got = 0.0_ep
            do src_rank = 0, my_nproc - 1
                allocate(buf(nb, nrhs))
                if (src_rank == my_rank) then
                    buf = B_full
                else
                    call mpi_recv(buf, nb * nrhs * bytes_per_elem, mpi_byte, &
                                  src_rank, tag, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    block
                        integer :: blk, owner
                        blk = (ig - 1) / nb; owner = blk
                        if (owner == src_rank) then
                            jl = mod(ig - 1, nb) + 1
                            B_full_got(ig, 1:nrhs) = buf(jl, 1:nrhs)
                        end if
                    end block
                end do
                if (src_rank == my_rank) then
                    buf = B_half
                else
                    call mpi_recv(buf, nb * nrhs * bytes_per_elem, mpi_byte, &
                                  src_rank, tag + 1, mpi_comm_world, &
                                  mpi_status_ignore, ierr)
                end if
                do ig = 1, n
                    block
                        integer :: blk, owner
                        blk = (ig - 1) / nb; owner = blk
                        if (owner == src_rank) then
                            jl = mod(ig - 1, nb) + 1
                            B_half_got(ig, 1:nrhs) = buf(jl, 1:nrhs)
                        end if
                    end block
                end do
                deallocate(buf)
            end do
        else
            peer = 0
            call mpi_send(B_full, nb * nrhs * bytes_per_elem, mpi_byte, &
                          peer, tag, mpi_comm_world, ierr)
            call mpi_send(B_half, nb * nrhs * bytes_per_elem, mpi_byte, &
                          peer, tag + 1, mpi_comm_world, ierr)
        end if

        if (my_rank == 0) then
            err = max_rel_err_mat(B_half_got, B_full_got)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_full_got, B_half_got)
        end if
        deallocate(d_loc, e_loc, B_full, B_half, work, af, &
                   d_glob, e_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpttrsv
