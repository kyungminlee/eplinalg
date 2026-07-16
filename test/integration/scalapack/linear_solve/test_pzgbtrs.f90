program test_pzgbtrs
    ! Complex general banded factor (pzgbtrf) + solve (pzgbtrs) split.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgesv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_nproc, &
                                my1d_context, my1d_npcol, my1d_col, &
                                descinit_1d
    use target_scalapack, only: target_name, target_eps, &
                                target_pzgbtrf, target_pzgbtrs
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: bwl = 2, bwu = 2, lda = 2*bwl + 2*bwu + 1
    integer :: i, n, nb, info, info_ref, lwork, laf, ig, jg, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer, parameter :: tag = 7751
    integer :: desca(9), descb(9)
    complex(ep), allocatable :: A_full(:,:), B_glob(:,:)
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: af(:)
    complex(ep), allocatable :: B_got(:,:), B_ref(:,:), buf(:,:)
    complex(ep), allocatable :: work(:)
    integer,     allocatable :: ipiv_loc(:), ipiv_ref(:)
    real(ep) :: err, tol
    complex(ep) :: wopt(1)
    character(len=48) :: label

    call grid_init()
    call report_init('pzgbtrs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        nb = (n + my_nproc - 1) / my_nproc
        allocate(A_full(n, n), B_glob(n, nrhs))
        A_full = (0.0_ep, 0.0_ep)
        block
            integer :: k, j, s
            do j = 1, n
                do k = max(1, j - bwu), min(n, j + bwl)
                    s = 20101 + 31 * i + 7 * j + 11 * k
                    A_full(k, j) = cmplx(real(mod(s, 1009), ep) / 1009.0_ep, &
                                         real(mod(s + 17, 1013), ep) / 1013.0_ep, ep)
                end do
                A_full(j, j) = A_full(j, j) + cmplx(real(2 * (bwl + bwu + 1), ep), 0.0_ep, ep)
                do k = 1, nrhs
                    s = 20201 + 31 * i + 7 * j + 11 * k
                    B_glob(j, k) = cmplx(real(mod(s, 1019), ep) / 1019.0_ep, &
                                         real(mod(s + 17, 1021), ep) / 1021.0_ep, ep)
                end do
            end do
        end block

        allocate(A_loc(lda, nb))
        A_loc = (0.0_ep, 0.0_ep)
        do jg = 1, n
            block
                integer :: blk, owner, ig2
                blk = (jg - 1) / nb; owner = blk
                if (owner == my1d_col) then
                    jl = mod(jg - 1, nb) + 1
                    do ig2 = max(1, jg - bwu), min(n, jg + bwl)
                        A_loc(bwl + 2 * bwu + 1 + ig2 - jg, jl) = A_full(ig2, jg)
                    end do
                end if
            end block
        end do

        lldB = nb
        allocate(B_loc(lldB, nrhs))
        B_loc = (0.0_ep, 0.0_ep)
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

        laf = (nb + bwu) * (bwl + bwu) + 6 * (bwl + bwu) * (bwl + 2 * bwu) + 16
        allocate(af(laf))
        ! PZGBTRF docs say IPIV size >= DESCA(NB), but the divide-and-
        ! conquer merge in pzgbtrf.f:1045 writes IPIV(LN+1 .. LN+BM+BMN)
        ! with LN <= NB-BW and BM+BMN up to 2*(BW+BWU); pad to cover it.
        allocate(ipiv_loc(nb + 2 * (bwl + bwu)))

        call target_pzgbtrf(n, bwl, bwu, A_loc, 1, desca, ipiv_loc, &
                            af, laf, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call target_pzgbtrf(n, bwl, bwu, A_loc, 1, desca, ipiv_loc, &
                            af, laf, work, lwork, info)
        deallocate(work)

        call target_pzgbtrs('N', n, bwl, bwu, nrhs, A_loc, 1, desca, ipiv_loc, &
                            B_loc, 1, descb, af, laf, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call target_pzgbtrs('N', n, bwl, bwu, nrhs, A_loc, 1, desca, ipiv_loc, &
                            B_loc, 1, descb, af, laf, work, lwork, info)

        bytes_per_elem = storage_size((0.0_ep, 0.0_ep)) / 8
        if (my_rank == 0) then
            allocate(B_got(n, nrhs))
            B_got = (0.0_ep, 0.0_ep)
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
            allocate(B_ref(n, nrhs), ipiv_ref(n))
            B_ref = B_glob
            call zgesv(n, nrhs, A_full, n, ipiv_ref, B_ref, n, info_ref)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0,a,i0,a,i0)') &
                'n=', n, ',bwl=', bwl, ',bwu=', bwu, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_ref, B_got, ipiv_ref)
        end if
        deallocate(A_loc, B_loc, work, af, ipiv_loc, A_full, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgbtrs
