program test_pddbtrsv
    ! Validate diagonally-dominant banded half-solves (real LU) by
    ! composing them to recover the full pddbtrs result.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_col, descinit_1d
    use mpi
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pddbtrf, target_pddbtrs, target_pddbtrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: bwl = 2, bwu = 2, lda = bwl + bwu + 1
    integer, parameter :: tag = 8861
    integer :: i, n, nb, info, lwork, laf, ig, jg, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_full_dense(:,:), B_glob(:,:)
    real(ep), allocatable :: A_loc(:,:), B_full(:,:), B_half(:,:)
    real(ep), allocatable :: af(:), af_save(:), work(:)
    real(ep), allocatable :: B_full_got(:,:), B_half_got(:,:), buf(:,:)
    real(ep) :: err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pddbtrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc
        allocate(A_full_dense(n, n), B_glob(n, nrhs))
        A_full_dense = 0.0_ep
        block
            integer :: k, j, s
            do j = 1, n
                do k = max(1, j - bwu), min(n, j + bwl)
                    s = 33601 + 31 * i + 7 * j + 11 * k
                    A_full_dense(k, j) = real(mod(s, 1009), ep) / 1009.0_ep
                end do
                A_full_dense(j, j) = A_full_dense(j, j) + real(2 * (bwl + bwu + 1), ep)
                do k = 1, nrhs
                    s = 33611 + 31 * i + 7 * j + 11 * k
                    B_glob(j, k) = real(mod(s, 1019), ep) / 1019.0_ep
                end do
            end do
        end block

        allocate(A_loc(lda, nb)); A_loc = 0.0_ep
        do jg = 1, n
            block
                integer :: blk, owner, ig2
                blk = (jg - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(jg - 1, nb) + 1
                    do ig2 = max(1, jg - bwu), min(n, jg + bwl)
                        A_loc(bwu + 1 + ig2 - jg, jl) = A_full_dense(ig2, jg)
                    end do
                end if
            end block
        end do

        lldB = nb
        allocate(B_full(lldB, nrhs), B_half(lldB, nrhs))
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

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, lda,  info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        laf = nb * (bwl + bwu) + 6 * max(bwl, bwu) * max(bwl, bwu) + 16
        allocate(af(laf), af_save(laf))

        call target_pddbtrf(n, bwl, bwu, A_loc, 1, desca, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pddbtrf(n, bwl, bwu, A_loc, 1, desca, af, laf, work, lwork, info)
        deallocate(work)

        ! Snapshot AF after the factor; both the full-trs path and the
        ! half-trsv path read it as unmodified factor data, but at least
        ! the trs internals can perturb internal scratch in AF on certain
        ! ranks. Restore between the two paths.
        af_save = af

        ! Half-solve composition first (trsv(L,N) then trsv(U,N) ≡ trs).
        ! Allocate the verified upper-bound workspace up front (the TODO
        ! note explains why upstream's pqdbtrsv WORK_SIZE_MIN under-allocates).
        lwork = nrhs * max(bwl, bwu) + (nrhs - 1) * (bwl + bwu)
        lwork = max(1, lwork)
        allocate(work(lwork))
        call target_pddbtrsv('L', 'N', n, bwl, bwu, nrhs, A_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)
        call target_pddbtrsv('U', 'N', n, bwl, bwu, nrhs, A_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        ! Now the full-solve oracle on the snapshotted AF.
        af = af_save
        call target_pddbtrs('N', n, bwl, bwu, nrhs, A_loc, 1, desca, &
                            B_full, 1, descb, af, laf, wopt, -1, info)
        ! Upstream pddbtrs LWMIN under-reports (see TODO entry on
        ! PDDBTRS / PZDBTRS workspace under-allocation); take the
        ! verified upper bound used in test_pddbtrs.f90.
        lwork = max(1, int(wopt(1)), &
                    nrhs * max(bwl, bwu) + (nrhs - 1) * (bwl + bwu))
        allocate(work(lwork))
        call target_pddbtrs('N', n, bwl, bwu, nrhs, A_loc, 1, desca, &
                            B_full, 1, descb, af, laf, work, lwork, info)

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
            write(label, '(a,i0,a,i0,a,i0,a,i0)') &
                'n=', n, ',bwl=', bwl, ',bwu=', bwu, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_full_got, B_half_got)
        end if
        deallocate(A_loc, B_full, B_half, work, af, af_save, A_full_dense, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pddbtrsv
