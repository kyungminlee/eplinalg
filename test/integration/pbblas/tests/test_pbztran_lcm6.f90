! Test PBZTRAN — non-square 3×2 grid exercising MIN(LCMP,CEIL(N,NB)) > 1
! in the ADIST='R' path, with TRANS='C' (conjugate transpose).
!
! Closes the "2×3 / 3×2 grids with MIN(LCMQ,CEIL(M,NB)) > 1" gap in
! tests/pbblas/TODO.md (mirror to test_pbdtran_lcm6, which covers the
! ADIST='C' / 2×3 half).
!
! On a 3×2 grid (NPROW=3, NPCOL=2): LCM(3,2)=6, LCMP=2, LCMQ=3.
! With N=12, NB=4 → CEIL(N,NB)=3 → MIN(LCMP=2, 3)=2, which exercises
! the LCMP-cycle DO 30 I = 0, LCMP-1 dispatch and the LCMQ.GT.1
! PBDTRGET / PBDTRSRT path in pbztran.f.
!
! Requires connected MPI (6 ranks). Skips cleanly when grid_init_shape
! detects a sandbox-style unconnected world.
program test_pbztran_lcm6
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_mat_z
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_matrix_complex
    use pblas_grid,         only: grid_init_shape, grid_exit, my_rank, &
                                  my_context, my_nprow, my_npcol, &
                                  my_row, my_col, numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbztran
    use mpi
    implicit none

    integer, parameter :: nb = 4
    integer, parameter :: m  = nb
    integer, parameter :: n  = 12
    integer, parameter :: iarow = 0, iacol = 0, icrow = 0, iccol = 0
    integer :: nq_loc, np_loc, lda, ldc, lwork
    integer :: bytes_per_elem, ig, jg, src_rank, peer, locm_C, locn_C, ierr2, pr
    integer, parameter :: tag_C = 7561
    complex(ep), allocatable :: A_glob(:,:), A_loc(:,:), C_loc(:,:), C0_loc(:,:)
    complex(ep), allocatable :: C_glob_got(:,:), C_glob_ref(:,:), C0_glob(:,:)
    complex(ep), allocatable :: buf2d(:,:)
    complex(ep) :: beta
    real(ep) :: err, tol
    character(len=96) :: label

    ! 3×2 grid: NPROW=3, NPCOL=2 ⇒ LCMP=2, LCMQ=3.
    call grid_init_shape(3, 2)
    call report_init('pbztran_lcm6', target_name, my_rank)

    if (my_context < 0) then
        if (my_rank == 0) then
            call report_case('skipped: nproc /= 6 (need connected MPI 6-rank world)', &
                             0.0_ep, 0.0_ep)
        end if
        call report_finalize()
        call grid_exit()
        return
    end if

    nq_loc = numroc_local(n, nb, my_col, iacol, my_npcol)
    np_loc = numroc_local(n, nb, my_row, icrow, my_nprow)
    lda = max(1, m)
    ldc = max(1, np_loc)
    lwork = max(64, 4 * n * m)

    call gen_matrix_complex(m, n, A_glob, seed = 4921)
    allocate(A_loc(lda, max(1, nq_loc)))
    A_loc = (0.0_ep, 0.0_ep)
    if (my_row == iarow .and. nq_loc > 0) then
        do jg = 1, n
            block
                integer :: zero_based, blk, owner_c, jl
                zero_based = jg - 1
                blk = zero_based / nb
                owner_c = mod(blk + iacol, my_npcol)
                if (owner_c == my_col) then
                    jl = (blk / my_npcol) * nb + mod(zero_based, nb) + 1
                    A_loc(1:m, jl) = A_glob(1:m, jg)
                end if
            end block
        end do
    end if

    call gen_matrix_complex(n, m, C0_glob, seed = 4931)
    allocate(C_loc(ldc, max(1, m)), C0_loc(ldc, max(1, m)))
    C_loc  = (0.0_ep, 0.0_ep)
    C0_loc = (0.0_ep, 0.0_ep)
    if (my_col == iccol .and. np_loc > 0) then
        do ig = 1, n
            block
                integer :: zero_based, blk, owner_r, il
                zero_based = ig - 1
                blk = zero_based / nb
                owner_r = mod(blk + icrow, my_nprow)
                if (owner_r == my_row) then
                    il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                    C_loc (il, 1:m) = C0_glob(ig, 1:m)
                    C0_loc(il, 1:m) = C0_glob(ig, 1:m)
                end if
            end block
        end do
    end if

    beta = (0.25_ep, 0.1_ep)
    call target_pbztran(my_context, 'R', 'C', m, n, nb, A_loc, lda, beta, &
                        C_loc, ldc, iarow, iacol, icrow, iccol, &
                        ncolsA = max(1, nq_loc), ncolsC = m, lwork = lwork)

    bytes_per_elem = storage_size((0.0_ep, 0.0_ep)) / 8
    if (my_rank == 0) then
        allocate(C_glob_got(n, m))
        C_glob_got = (0.0_ep, 0.0_ep)
        do pr = 0, my_nprow - 1
            locm_C = numroc_local(n, nb, pr, icrow, my_nprow)
            locn_C = m
            if (locm_C == 0) cycle
            allocate(buf2d(locm_C, locn_C))
            if (pr == my_row .and. my_col == iccol) then
                buf2d = C_loc(1:locm_C, 1:locn_C)
            else
                src_rank = pr * my_npcol + iccol
                call mpi_recv(buf2d, locm_C * locn_C * bytes_per_elem, &
                              mpi_byte, src_rank, tag_C, mpi_comm_world, &
                              mpi_status_ignore, ierr2)
            end if
            do ig = 1, n
                block
                    integer :: zero_based, blk, owner_r, il
                    zero_based = ig - 1
                    blk = zero_based / nb
                    owner_r = mod(blk + icrow, my_nprow)
                    if (owner_r == pr) then
                        il = (blk / my_nprow) * nb + mod(zero_based, nb) + 1
                        C_glob_got(ig, 1:m) = buf2d(il, 1:m)
                    end if
                end block
            end do
            deallocate(buf2d)
        end do
    else
        if (my_col == iccol .and. np_loc > 0) then
            peer = 0
            call mpi_send(C_loc, np_loc * m * bytes_per_elem, mpi_byte, peer, &
                          tag_C, mpi_comm_world, ierr2)
        end if
    end if

    if (my_rank == 0) then
        allocate(C_glob_ref(n, m))
        C_glob_ref = transpose(conjg(A_glob)) + beta * C0_glob
        err = max_rel_err_mat_z(C_glob_got, C_glob_ref)
        tol = 64.0_ep * real(m * n, ep) * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0,a,i0)') &
            'grid=', my_nprow, 'x', my_npcol, &
            ',adist=R,LCMP=2,m=', m, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(C_glob_got, C_glob_ref)
    end if

    deallocate(A_glob, A_loc, C_loc, C0_loc, C0_glob)

    call report_finalize()
    call grid_exit()
end program test_pbztran_lcm6
