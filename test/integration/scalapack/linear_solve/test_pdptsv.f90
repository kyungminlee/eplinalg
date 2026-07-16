program test_pdptsv
    ! ScaLAPACK 1D-distributed positive-definite tridiagonal solve.
    ! Mirrors test_pddtsv but uses the symmetric form: D is the
    ! diagonal and E is the (single) off-diagonal.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dptsv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_nproc, &
                                my1d_context, my1d_npcol, my1d_col, &
                                descinit_1d
    use test_data,        only: gen_vector_quad, gen_matrix_quad
    use target_scalapack, only: target_name, target_eps, target_pdptsv
    use mpi
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer :: i, n, nb, info, info_ref, lwork, ig, jl, lldB
    integer :: bytes_per_elem, ierr, src_rank, peer
    integer, parameter :: tag = 7621
    integer :: desca(9), descb(9)
    real(ep), allocatable :: d_glob(:), e_glob(:), B_glob(:,:)
    real(ep), allocatable :: d_loc(:),  e_loc(:),  B_loc(:,:)
    real(ep), allocatable :: B_got(:,:), B_ref(:,:), buf(:,:)
    real(ep), allocatable :: d_ref(:), e_ref(:), work(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdptsv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        nb = (n + my_nproc - 1) / my_nproc

        call gen_vector_quad(n,     d_glob, seed = 18101 + 31*i)
        call gen_vector_quad(n,     e_glob, seed = 18111 + 31*i)
        call gen_matrix_quad(n, nrhs, B_glob, seed = 18121 + 31*i)
        ! Make D positive and dominant; clamp E to keep PD.
        d_glob = abs(d_glob) + real(2 * n, ep)
        e_glob = e_glob * 0.1_ep   ! small off-diagonal
        e_glob(n) = 0.0_ep         ! E(n) ignored / unused

        allocate(d_loc(nb), e_loc(nb))
        d_loc = 0.0_ep; e_loc = 0.0_ep
        do ig = 1, n
            block
                integer :: blk, owner
                blk = (ig - 1) / nb
                owner = blk
                if (owner == my1d_col) then
                    jl = mod(ig - 1, nb) + 1
                    d_loc(jl) = d_glob(ig)
                    e_loc(jl) = e_glob(ig)
                end if
            end block
        end do

        lldB = nb
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

        call descinit_1d(desca, 501, n, nb, 0, my1d_context, nb,   info)
        call descinit_1d(descb, 502, n, nb, 0, my1d_context, lldB, info)

        ! pdptsv LWORK >= 12*NPCOL + 3*NB + max((10+2*min(100,NRHS))*NPCOL+4*NRHS, 8*NPCOL).
        lwork = 12 * my1d_npcol + 3 * nb &
              + max((10 + 2 * min(100, nrhs)) * my1d_npcol + 4 * nrhs, 8 * my1d_npcol)
        allocate(work(max(1, lwork)))

        call target_pdptsv(n, nrhs, d_loc, e_loc, 1, desca, &
                           B_loc, 1, descb, work, lwork, info)

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
            allocate(d_ref(n), e_ref(n), B_ref(n, nrhs))
            ! dptsv expects e(1..n-1) for the off-diagonal; pdptsv
            ! also uses e(1..n) with e(n) unused. They line up.
            d_ref = d_glob
            e_ref(1:n-1) = e_glob(1:n-1)
            e_ref(n) = 0.0_ep
            B_ref = B_glob
            call dptsv(n, nrhs, d_ref, e_ref, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(d_ref, e_ref, B_ref, B_got)
        end if
        deallocate(d_loc, e_loc, B_loc, work, d_glob, e_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdptsv
