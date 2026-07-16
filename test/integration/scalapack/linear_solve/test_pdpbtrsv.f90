program test_pdpbtrsv
    ! Validate SPD banded half-solves (UPLO='U' factor A = U^T*U) by
    ! composing them to recover the full pdpbtrs result. A^-1 * B factors
    ! as U^-1 * (U^-T * B), so apply pdpbtrsv('U','T') first then
    ! pdpbtrsv('U','N').
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_nproc, &
                                 my1d_context, my1d_col, descinit_1d
    use mpi
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pdpbtrf, target_pdpbtrs, target_pdpbtrsv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: bw = 2, lda = bw + 1
    integer, parameter :: tag = 8841
    integer :: i, n, nb, info, lwork, laf, ig, jg, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_full_dense(:,:), M(:,:), B_glob(:,:)
    real(ep), allocatable :: A_loc(:,:), B_full(:,:), B_half(:,:)
    real(ep), allocatable :: af(:), work(:)
    real(ep), allocatable :: B_full_got(:,:), B_half_got(:,:), buf(:,:)
    real(ep) :: err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pdpbtrsv', target_name, my_rank)

    do i = 1, size(ns)
        n  = ns(i)
        nb = (n + my_nproc - 1) / my_nproc
        allocate(A_full_dense(n, n), M(n, n), B_glob(n, nrhs))
        block
            integer :: k, j, s
            M = 0.0_ep
            do j = 1, n
                do k = max(1, j - bw), min(n, j + bw)
                    s = 33501 + 31 * i + 7 * j + 11 * k
                    M(k, j) = real(mod(s, 1009), ep) / 1009.0_ep
                end do
                do k = 1, nrhs
                    s = 33511 + 31 * i + 7 * j + 11 * k
                    B_glob(j, k) = real(mod(s, 1019), ep) / 1019.0_ep
                end do
            end do
        end block
        A_full_dense = matmul(transpose(M), M)
        block
            integer :: k, j
            do k = 1, n
                A_full_dense(k, k) = A_full_dense(k, k) + real(n, ep)
            end do
            do j = 1, n
                do k = 1, n
                    if (abs(k - j) > bw) A_full_dense(k, j) = 0.0_ep
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
                    do ig2 = max(1, jg - bw), jg
                        A_loc(bw + 1 + ig2 - jg, jl) = A_full_dense(ig2, jg)
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

        laf = (nb + 2 * bw) * bw + 16
        allocate(af(laf))

        call target_pdpbtrf('U', n, bw, A_loc, 1, desca, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpbtrf('U', n, bw, A_loc, 1, desca, af, laf, work, lwork, info)
        deallocate(work)

        call target_pdpbtrs('U', n, bw, nrhs, A_loc, 1, desca, &
                            B_full, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpbtrs('U', n, bw, nrhs, A_loc, 1, desca, &
                            B_full, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        ! Apply U^-T then U^-1 to B_half.
        call target_pdpbtrsv('U', 'T', n, bw, nrhs, A_loc, 1, desca, &
                             B_half, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpbtrsv('U', 'T', n, bw, nrhs, A_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, lwork, info)
        deallocate(work)

        allocate(work(1))
        call target_pdpbtrsv('U', 'N', n, bw, nrhs, A_loc, 1, desca, &
                             B_half, 1, descb, af, laf, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdpbtrsv('U', 'N', n, bw, nrhs, A_loc, 1, desca, &
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
            write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',bw=', bw, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_full_got, B_half_got)
        end if
        deallocate(A_loc, B_full, B_half, work, af, A_full_dense, M, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpbtrsv
