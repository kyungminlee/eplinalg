! Test PBDTRAN — non-square 2×3 grid exercising MIN(LCMQ,CEIL(M,NB)) > 1
! in the ADIST='C' path.
!
! Closes the "2×3 / 3×2 grids with MIN(LCMQ,CEIL(M,NB)) > 1" gap in
! tests/pbblas/TODO.md. The companion test_pbztran_lcm6 covers the
! mirror ADIST='R' path on a 3×2 grid.
!
! On a 2×3 grid (NPROW=2, NPCOL=3): LCM(2,3)=6, LCMP=3, LCMQ=2.
! With M=12, NB=4 → CEIL(M,NB)=3 → MIN(LCMQ=2, 3)=2, which activates
! the LCM-cycle dispatch loop (DO 20 I = 0, MIN(LCM,CEIL(M,NB))-1) and
! the LCMQ.GT.1 / PBDTRGET / PBDTRSRT branches inside pbdtran.f that
! a 2×2 grid (LCMP=LCMQ=1) collapses past.
!
! Requires connected MPI (6 ranks). Like test_pbdtrnv_lcm, this skips
! cleanly when grid_init_shape sees a sandbox-style unconnected world.
! Use the linux-impi preset (PBBLAS_TEST_NPROC_LCM6=6) to actually
! exercise the path.
program test_pbdtran_lcm6
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_mat
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_matrix_quad
    use pblas_grid,         only: grid_init_shape, grid_exit, my_rank, &
                                  my_context, my_nprow, my_npcol, &
                                  my_row, my_col, numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbdtran
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: m  = 12
    integer, parameter :: n  = nb
    integer, parameter :: iarow = 0, iacol = 0, icrow = 0, iccol = 0
    integer :: mp, mq, lda, ldc, lwork
    integer :: bytes_per_elem, ig, jg, src_rank, peer, locn_C, ierr2, pc
    integer, parameter :: tag_C = 7551
    real(ep), allocatable :: A_glob(:,:), A_loc(:,:), C_loc(:,:), C0_loc(:,:)
    real(ep), allocatable :: C_glob_got(:,:), C_glob_ref(:,:), C0_glob(:,:)
    real(ep), allocatable :: buf2d(:,:)
    real(ep) :: beta, err, tol
    character(len=96) :: label

    ! 2×3 grid: NPROW=2, NPCOL=3 ⇒ LCMP=3, LCMQ=2.
    call grid_init_shape(2, 3)
    call report_init('pbdtran_lcm6', target_name, my_rank)

    if (my_context < 0) then
        if (my_rank == 0) then
            call report_case('skipped: nproc /= 6 (need connected MPI 6-rank world)', &
                             0.0_ep, 0.0_ep)
        end if
        call report_finalize()
        call grid_exit()
        return
    end if

    mp = numroc_local(m, nb, my_row, iarow, my_nprow)
    mq = numroc_local(m, nb, my_col, iccol, my_npcol)
    lda = max(1, mp)
    ldc = max(1, n)
    lwork = max(64, 4 * n * m)

    call gen_matrix_quad(m, n, A_glob, seed = 4901)
    allocate(A_loc(lda, max(1, n)))
    A_loc = 0.0_ep
    if (my_col == iacol) then
        do ig = 1, m
            block
                integer :: zero_based, blk, owner_r, il
                zero_based = ig - 1
                blk = zero_based / nb
                owner_r = mod(blk + iarow, my_nprow)
                if (owner_r == my_row) then
                    il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                    A_loc(il, 1:n) = A_glob(ig, 1:n)
                end if
            end block
        end do
    end if

    call gen_matrix_quad(n, m, C0_glob, seed = 4911)
    allocate(C_loc(ldc, max(1, mq)), C0_loc(ldc, max(1, mq)))
    C_loc = 0.0_ep
    C0_loc = 0.0_ep
    if (my_row == icrow .and. mq > 0) then
        do jg = 1, m
            block
                integer :: zero_based, blk, owner_c, jl
                zero_based = jg - 1
                blk = zero_based / nb
                owner_c = mod(blk + iccol, my_npcol)
                if (owner_c == my_col) then
                    jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                    C_loc (1:n, jl) = C0_glob(1:n, jg)
                    C0_loc(1:n, jl) = C0_glob(1:n, jg)
                end if
            end block
        end do
    end if

    beta = 0.5_ep
    call target_pbdtran(my_context, 'C', 'T', m, n, nb, A_loc, lda, beta, &
                        C_loc, ldc, iarow, iacol, icrow, iccol, &
                        ncolsA = n, ncolsC = max(1, mq), lwork = lwork)

    bytes_per_elem = storage_size(0.0_ep) / 8
    if (my_rank == 0) then
        allocate(C_glob_got(n, m))
        C_glob_got = 0.0_ep
        do pc = 0, my_npcol - 1
            locn_C = numroc_local(m, nb, pc, iccol, my_npcol)
            if (locn_C == 0) cycle
            allocate(buf2d(n, locn_C))
            if (pc == my_col .and. my_row == icrow) then
                buf2d = C_loc(1:n, 1:locn_C)
            else
                src_rank = icrow * my_npcol + pc
                call mpi_recv(buf2d, n * locn_C * bytes_per_elem, &
                              mpi_byte, src_rank, tag_C, mpi_comm_world, &
                              mpi_status_ignore, ierr2)
            end if
            do jg = 1, m
                block
                    integer :: zero_based, blk, owner_c, jl
                    zero_based = jg - 1
                    blk = zero_based / nb
                    owner_c = mod(blk + iccol, my_npcol)
                    if (owner_c == pc) then
                        jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                        C_glob_got(1:n, jg) = buf2d(1:n, jl)
                    end if
                end block
            end do
            deallocate(buf2d)
        end do
    else
        if (my_row == icrow .and. mq > 0) then
            peer = 0
            call mpi_send(C_loc, n * mq * bytes_per_elem, mpi_byte, peer, &
                          tag_C, mpi_comm_world, ierr2)
        end if
    end if

    if (my_rank == 0) then
        allocate(C_glob_ref(n, m))
        C_glob_ref = transpose(A_glob) + beta * C0_glob
        err = max_rel_err_mat(C_glob_got, C_glob_ref)
        tol = 64.0_ep * real(m * n, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') &
            'grid=', my_nprow, 'x', my_npcol, &
            ',adist=C,LCMQ=2,m=', m, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(C_glob_got, C_glob_ref)
    end if

    deallocate(A_glob, A_loc, C_loc, C0_loc, C0_glob)

    call report_finalize()
    call grid_exit()
end program test_pbdtran_lcm6
