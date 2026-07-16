program test_pdpbtrs
    ! Symmetric-PD banded factor (pdpbtrf) + solve (pdpbtrs) split.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dposv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_nproc, &
                                my1d_context, my1d_npcol, my1d_col, &
                                descinit_1d
    use target_scalapack, only: target_name, target_eps, &
                                target_pdpbtrf, target_pdpbtrs
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: bw = 2, lda = bw + 1
    integer :: i, n, nb, info, info_ref, lwork, laf, ig, jg, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer, parameter :: tag = 7781
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_full(:,:), B_glob(:,:), M(:,:)
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    real(ep), allocatable :: af(:)
    real(ep), allocatable :: B_got(:,:), B_ref(:,:), buf(:,:), A_ref(:,:)
    real(ep), allocatable :: work(:)
    real(ep) :: err, tol, wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pdpbtrs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        nb = (n + my_nproc - 1) / my_nproc
        allocate(A_full(n, n), B_glob(n, nrhs), M(n, n))
        block
            integer :: k, j, s
            do j = 1, n
                do k = 1, n
                    M(k, j) = 0.0_ep
                end do
                do k = max(1, j - bw), min(n, j + bw)
                    s = 20701 + 31 * i + 7 * j + 11 * k
                    M(k, j) = real(mod(s, 1009), ep) / 1009.0_ep
                end do
                do k = 1, nrhs
                    s = 20801 + 31 * i + 7 * j + 11 * k
                    B_glob(j, k) = real(mod(s, 1019), ep) / 1019.0_ep
                end do
            end do
        end block
        A_full = matmul(transpose(M), M)
        block
            integer :: k, j
            do k = 1, n
                A_full(k, k) = A_full(k, k) + real(n, ep)
            end do
            do j = 1, n
                do k = 1, n
                    if (abs(k - j) > bw) A_full(k, j) = 0.0_ep
                end do
            end do
        end block

        allocate(A_loc(lda, nb))
        A_loc = 0.0_ep
        do jg = 1, n
            block
                integer :: blk, owner, ig2
                blk = (jg - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(jg - 1, nb) + 1
                    do ig2 = max(1, jg - bw), jg
                        A_loc(bw + 1 + ig2 - jg, jl) = A_full(ig2, jg)
                    end do
                end if
            end block
        end do

        lldB = nb
        allocate(B_loc(lldB, nrhs))
        B_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    B_loc(jl, 1:nrhs) = B_glob(ig, 1:nrhs)
                end if
            end block
        end do

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, lda,  info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        laf = (nb + 2 * bw) * bw + 16
        allocate(af(laf))

        call target_pdpbtrf('U', n, bw, A_loc, 1, desca, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpbtrf('U', n, bw, A_loc, 1, desca, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pdpbtrs('U', n, bw, nrhs, A_loc, 1, desca, &
                            B_loc, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdpbtrs('U', n, bw, nrhs, A_loc, 1, desca, &
                            B_loc, 1, descb, af, laf, work, lwork, info)

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
                        blk = (ig - 1) / nb; owner = blk
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
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_full
            block
                integer :: k, j
                do j = 1, n
                    do k = 1, j - 1
                        A_ref(j, k) = A_ref(k, j)
                    end do
                end do
            end block
            B_ref = B_glob
            call dposv('U', n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',bw=', bw, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, B_got)
        end if
        deallocate(A_loc, B_loc, work, af, A_full, B_glob, M)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpbtrs
