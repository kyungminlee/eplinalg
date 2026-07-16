program test_pdtrcon
    ! 1-norm condition number estimate for a real triangular matrix.
    ! Same Hager-Higham slop as pdgecon/pdpocon — compare ScaLAPACK
    ! rcond against the LAPACK serial estimate in orders of magnitude
    ! (tolerance = 2 orders), not in target_eps.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dtrcon
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdtrcon
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_ref(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: iwork(:), iwork_ref(:)
    real(ep) :: rcond, rcond_ref, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdtrcon', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 31101 + 31*i)

        ! Force upper-triangular and diagonally biased.
        do jg = 1, n
            do ig = jg + 1, n
                A_glob(ig, jg) = 0.0_ep
            end do
            A_glob(jg, jg) = A_glob(jg, jg) + real(2 * n, ep)
        end do
        if (size(A_loc, 1) > 0 .and. size(A_loc, 2) > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_glob(ig, jg)
                end do
            end do
        end if

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(work(1), iwork(1))
        call target_pdtrcon('1', 'U', 'N', n, A_loc, 1, 1, desca, rcond, &
                            work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdtrcon('1', 'U', 'N', n, A_loc, 1, 1, desca, rcond, &
                            work, lwork, iwork, liwork, info)

        if (my_rank == 0) then
            allocate(A_ref(n, n), work_ref(3 * n), iwork_ref(n))
            A_ref = A_glob
            call dtrcon('1', 'U', 'N', n, A_ref, n, rcond_ref, &
                        work_ref, iwork_ref, info_ref)
            if (rcond /= 0.0_ep .and. rcond_ref /= 0.0_ep) then
                err = abs(log(rcond / rcond_ref)) / log(10.0_ep)
            else
                err = 0.0_ep
            end if
            tol = 2.0_ep
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, work_ref, iwork_ref)
        end if
        deallocate(A_loc, A_glob, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtrcon
