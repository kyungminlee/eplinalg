program test_pddttrsv
    ! Validate the general-tridiagonal half-solves by composing them to
    ! recover the full pddttrs result. PDDTTRF stores A = P*L*U with the
    ! permutation absorbed into the factor representation; PDDTTRSV
    ! halves with TRANS='N' apply L^-1 (UPLO='L') and U^-1 (UPLO='U')
    ! through the same factor.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_npcol, my1d_col, descinit_1d
    use test_data,         only: gen_vector_quad, gen_matrix_quad
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pddttrf, target_pddttrs, target_pddttrsv
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: tag = 8821
    integer :: i, n, nb, info, lwork, laf, ig, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer :: desca(9), descb(9)
    real(ep), allocatable :: dl_glob(:), d_glob(:), du_glob(:), B_glob(:,:)
    real(ep), allocatable :: dl_loc(:), d_loc(:), du_loc(:)
    real(ep), allocatable :: B_full(:,:), B_half(:,:)
    real(ep), allocatable :: af(:), work(:)
    real(ep), allocatable :: B_full_got(:,:), B_half_got(:,:), buf(:,:)
    real(ep) :: err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pddttrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n, dl_glob, seed = 33401 + 31*i)
        call gen_vector_quad(n, d_glob,  seed = 33411 + 31*i)
        call gen_vector_quad(n, du_glob, seed = 33421 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 33431 + 31*i)
        dl_glob = dl_glob * 0.05_ep
        du_glob = du_glob * 0.05_ep
        d_glob  = abs(d_glob) + real(2 * n, ep)
        dl_glob(1) = 0.0_ep
        du_glob(n) = 0.0_ep

        allocate(dl_loc(nb), d_loc(nb), du_loc(nb))
        dl_loc = 0.0_ep; d_loc = 0.0_ep; du_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    dl_loc(jl) = dl_glob(ig)
                    d_loc(jl)  = d_glob(ig)
                    du_loc(jl) = du_glob(ig)
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

        call target_pddttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddttrf(n, dl_loc, d_loc, du_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pddttrs('N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                            B_full, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddttrs('N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                            B_full, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        call target_pddttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_half, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddttrsv('L', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        allocate(work(1))
        call target_pddttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pddttrsv('U', 'N', n, nrhs, dl_loc, d_loc, du_loc, 1, desca, &
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
        deallocate(dl_loc, d_loc, du_loc, B_full, B_half, work, af, &
                   dl_glob, d_glob, du_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pddttrsv
