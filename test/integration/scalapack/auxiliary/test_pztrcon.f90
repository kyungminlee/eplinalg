program test_pztrcon
    ! 1-norm condition number estimate for a complex triangular matrix.
    ! Hager-Higham slop — compare orders of magnitude (tol = 2).
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: ztrcon
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z
    use target_scalapack,  only: target_name, target_eps, target_pztrcon
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, lrwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_ref(:,:)
    complex(ep), allocatable :: work(:), work_ref(:)
    real(ep), allocatable :: rwork(:), rwork_ref(:)
    real(ep) :: rcond, rcond_ref, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pztrcon', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 31201 + 31*i)

        do jg = 1, n
            do ig = jg + 1, n
                A_glob(ig, jg) = (0.0_ep, 0.0_ep)
            end do
            A_glob(jg, jg) = A_glob(jg, jg) + cmplx(real(2 * n, ep), 0.0_ep, ep)
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

        allocate(work(1), rwork(1))
        call target_pztrcon('1', 'U', 'N', n, A_loc, 1, 1, desca, rcond, &
                            work, -1, rwork, -1, info)
        lwork  = max(1, int(real(work(1))))
        lrwork = max(1, int(rwork(1)))
        deallocate(work, rwork)
        allocate(work(lwork), rwork(lrwork))
        call target_pztrcon('1', 'U', 'N', n, A_loc, 1, 1, desca, rcond, &
                            work, lwork, rwork, lrwork, info)

        if (my_rank == 0) then
            allocate(A_ref(n, n), work_ref(2 * n), rwork_ref(n))
            A_ref = A_glob
            call ztrcon('1', 'U', 'N', n, A_ref, n, rcond_ref, &
                        work_ref, rwork_ref, info_ref)
            if (rcond /= 0.0_ep .and. rcond_ref /= 0.0_ep) then
                err = abs(log(rcond / rcond_ref)) / log(10.0_ep)
            else
                err = 0.0_ep
            end if
            tol = 2.0_ep
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, work_ref, rwork_ref)
        end if
        deallocate(A_loc, A_glob, work, rwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrcon
