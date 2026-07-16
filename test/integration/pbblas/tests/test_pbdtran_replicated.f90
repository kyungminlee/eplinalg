! Test PBDTRAN — replicated-A paths.
!
! Closes the IACOL=-1 / IAROW=-1 coverage gap in tests/pbblas/TODO.md.
! pbdtran's replicated branches (line ~352 onward in pbdtran.f, the
! "When all column processors have a copy of the column block A" /
! mirror) are entirely separate from the IACOL>=0 / IAROW>=0 branches,
! using `PBDTRGET` / `PBDTRSRT` to gather diagonal blocks rather than
! the source-node-sends pattern.
!
! Default grid is 2×2: NPCOL=2 makes IACOL=-1 a real replication
! across two process columns (and ditto for IAROW=-1 over two rows).
! LCMP=LCMQ=1 on this grid so the inner LCM-cycle logic still
! collapses, but the *branches themselves* are exercised — that is
! the gap the TODO bullet referred to.
!
! Two cases:
!   ADIST='C', IACOL=-1, ICROW=0      — A replicated across columns
!   ADIST='R', IAROW=-1, ICCOL=0      — A replicated across rows
!
! C := A^T + beta*C0  (alpha is implicit; see test_pbdtran.f90).
program test_pbdtran_replicated
    use prec_kinds,         only: ep
    use compare,            only: max_rel_err_mat
    use pbblas_prec_report, only: report_init, report_case, report_finalize
    use test_data,          only: gen_matrix_quad
    use pblas_grid,         only: grid_init, grid_exit, my_rank, my_context, &
                                  my_nprow, my_npcol, my_row, my_col, &
                                  numroc_local
    use target_pbblas,      only: target_name, target_eps, target_pbdtran
    use mpi
    implicit none

    integer, parameter :: nb = 4

    call grid_init()
    call report_init('pbdtran_replicated', target_name, my_rank)

    call run_case_C_iacol_neg()
    call run_case_R_iarow_neg()

    call report_finalize()
    call grid_exit()

contains

    subroutine run_case_C_iacol_neg()
        ! ADIST='C', IACOL=-1: every process column holds A; A's rows
        ! are block-cyclic across NPROW (no IACOL routing).
        integer, parameter :: m = 8
        integer, parameter :: n = nb
        integer, parameter :: iarow = 0, iacol = -1, icrow = 0, iccol = 0
        integer :: mp, mq, lda, ldc, lwork
        integer :: bytes_per_elem, ig, jg, locm_C, locn_C, src_rank, peer
        integer, parameter :: tag_C = 7321
        real(ep), allocatable :: A_glob(:,:), A_loc(:,:), C_loc(:,:), C0_loc(:,:)
        real(ep), allocatable :: C_glob_got(:,:), C_glob_ref(:,:), C0_glob(:,:)
        real(ep), allocatable :: buf2d(:,:)
        real(ep) :: beta, err, tol
        character(len=80) :: label
        integer :: ierr2, pc

        mp = numroc_local(m, nb, my_row, iarow, my_nprow)
        mq = numroc_local(m, nb, my_col, iccol, my_npcol)
        lda = max(1, mp)
        ldc = max(1, n)
        ! Replicated-A WORK upper bound:
        !   N * CEIL(Mqb,LCMQ)*NB * MIN(LCMQ,CEIL(M,NB))
        ! N*M*MIN(NPROW,NPCOL) envelopes the formula on small grids.
        lwork = max(64, n * m * max(1, min(my_nprow, my_npcol)))

        call gen_matrix_quad(m, n, A_glob, seed = 4401)
        allocate(A_loc(lda, max(1, n)))
        A_loc = 0.0_ep
        ! IACOL=-1 ⇒ every column holds A.
        if (mp > 0) then
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

        call gen_matrix_quad(n, m, C0_glob, seed = 4501)
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
                locm_C = n
                locn_C = numroc_local(m, nb, pc, iccol, my_npcol)
                if (locn_C == 0) cycle
                allocate(buf2d(locm_C, locn_C))
                if (pc == my_col .and. my_row == icrow) then
                    buf2d = C_loc(1:locm_C, 1:locn_C)
                else
                    src_rank = icrow * my_npcol + pc
                    call mpi_recv(buf2d, locm_C * locn_C * bytes_per_elem, &
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
            write(label, '(a,a,i0,a,i0,a,i0)') 'adist=C,iacol=-1,', 'm=', m, &
                                                ',n=', n, ',nb=', nb
            call report_case(trim(label), err, tol)
            deallocate(C_glob_got, C_glob_ref)
        end if

        deallocate(A_glob, A_loc, C_loc, C0_loc, C0_glob)
    end subroutine run_case_C_iacol_neg

    subroutine run_case_R_iarow_neg()
        ! ADIST='R', IAROW=-1: every process row holds A; A's columns
        ! are block-cyclic across NPCOL.
        integer, parameter :: m = nb
        integer, parameter :: n = 8
        integer, parameter :: iarow = -1, iacol = 0, icrow = 0, iccol = 0
        integer :: nq_loc, np_loc, lda, ldc, lwork
        integer :: bytes_per_elem, ig, jg, locm_C, locn_C, src_rank, peer
        integer, parameter :: tag_C = 7331
        real(ep), allocatable :: A_glob(:,:), A_loc(:,:), C_loc(:,:), C0_loc(:,:)
        real(ep), allocatable :: C_glob_got(:,:), C_glob_ref(:,:), C0_glob(:,:)
        real(ep), allocatable :: buf2d(:,:)
        real(ep) :: beta, err, tol
        character(len=80) :: label
        integer :: ierr2, pr

        nq_loc = numroc_local(n, nb, my_col, iacol, my_npcol)
        np_loc = numroc_local(n, nb, my_row, icrow, my_nprow)
        lda = max(1, m)
        ldc = max(1, np_loc)
        lwork = max(64, n * m * max(1, min(my_nprow, my_npcol)))

        call gen_matrix_quad(m, n, A_glob, seed = 4601)
        allocate(A_loc(lda, max(1, nq_loc)))
        A_loc = 0.0_ep
        ! IAROW=-1 ⇒ every row holds A.
        if (nq_loc > 0) then
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

        call gen_matrix_quad(n, m, C0_glob, seed = 4701)
        allocate(C_loc(ldc, max(1, m)), C0_loc(ldc, max(1, m)))
        C_loc = 0.0_ep
        C0_loc = 0.0_ep
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

        beta = 0.25_ep
        call target_pbdtran(my_context, 'R', 'T', m, n, nb, A_loc, lda, beta, &
                            C_loc, ldc, iarow, iacol, icrow, iccol, &
                            ncolsA = max(1, nq_loc), ncolsC = m, lwork = lwork)

        bytes_per_elem = storage_size(0.0_ep) / 8
        if (my_rank == 0) then
            allocate(C_glob_got(n, m))
            C_glob_got = 0.0_ep
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
            C_glob_ref = transpose(A_glob) + beta * C0_glob
            err = max_rel_err_mat(C_glob_got, C_glob_ref)
            tol = 64.0_ep * real(m * n, ep) * target_eps
            write(label, '(a,a,i0,a,i0,a,i0)') 'adist=R,iarow=-1,', 'm=', m, &
                                                ',n=', n, ',nb=', nb
            call report_case(trim(label), err, tol)
            deallocate(C_glob_got, C_glob_ref)
        end if

        deallocate(A_glob, A_loc, C_loc, C0_loc, C0_glob)
    end subroutine run_case_R_iarow_neg

end program test_pbdtran_replicated
