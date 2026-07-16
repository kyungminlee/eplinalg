program test_pdrscl
    ! Vector reciprocal scaling: x := x / sa. Verify against the
    ! quad-computed analytic result.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdrscl
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info
    integer :: locm_x, lld_x
    integer :: descx(9)
    real(ep), allocatable :: x_loc(:,:), x_glob(:,:), x_got(:,:)
    real(ep) :: sa, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdrscl', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        ! Treat the vector as an N x 1 distributed matrix.
        call gen_distrib_matrix(n, 1, mb, nb, x_loc, x_glob, seed = 23201 + 31*i)
        sa = 7.0_ep + 0.1_ep * real(i, ep)

        locm_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locm_x)
        call descinit_local(descx, n, 1, mb, nb, 0, 0, my_context, lld_x, info)

        call target_pdrscl(n, sa, x_loc, 1, 1, descx, 1)
        call gather_matrix(n, 1, mb, nb, x_loc, x_got)

        if (my_rank == 0) then
            err = max_rel_err_vec(x_got(:, 1), x_glob(:, 1) / sa)
            tol = 32.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_got)
        end if
        deallocate(x_loc, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdrscl
