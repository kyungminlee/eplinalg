program test_pdswap
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dswap
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector, gather_vector
    use target_pblas,  only: target_name, target_eps, target_pdswap
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9), descy(9)
    real(ep), allocatable :: x_loc(:), y_loc(:)
    real(ep), allocatable :: x_glob(:), y_glob(:), x_got(:), y_got(:)
    real(ep), allocatable :: x_ref(:), y_ref(:)
    real(ep) :: err_x, err_y, err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pdswap', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector(n, mb, x_loc, x_glob, seed = 1701 + 7 * i)
        call gen_distrib_vector(n, mb, y_loc, y_glob, seed = 1801 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)
        call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld, info)

        call target_pdswap(n, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1)
        call gather_vector(n, mb, x_loc, x_got)
        call gather_vector(n, mb, y_loc, y_got)

        if (my_rank == 0) then
            allocate(x_ref(n), y_ref(n))
            x_ref = x_glob; y_ref = y_glob
            call dswap(n, x_ref, 1, y_ref, 1)
            err_x = max_rel_err_vec(x_got, x_ref)
            err_y = max_rel_err_vec(y_got, y_ref)
            err = max(err_x, err_y)
            ! Swap is exact — no arithmetic, just reads/writes.
            tol = 2.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, y_ref, x_got, y_got)
        end if
        deallocate(x_loc, y_loc, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdswap
