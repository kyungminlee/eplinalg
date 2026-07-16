program test_pzcopy
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zcopy
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z, gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzcopy
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9), descy(9)
    complex(ep), allocatable :: x_loc(:), y_loc(:)
    complex(ep), allocatable :: x_glob(:), y_glob(:), y_got(:), y_ref(:)
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzcopy', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2201 + 7 * i)
        call gen_distrib_vector_z(n, mb, y_loc, y_glob, seed = 2207 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)
        call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld, info)

        call target_pzcopy(n, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1)
        call gather_vector_z(n, mb, y_loc, y_got)

        if (my_rank == 0) then
            allocate(y_ref(n))
            call zcopy(n, x_glob, 1, y_ref, 1)
            err = max_rel_err_vec_z(y_got, y_ref)
            ! Copy is exact.
            tol = 2.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(y_ref, y_got)
        end if
        deallocate(x_loc, y_loc, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzcopy
