program test_pzdotc
    use prec_kinds,    only: ep
    use compare,       only: rel_err_scalar_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zdotc
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzdotc
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9), descy(9)
    complex(ep), allocatable :: x_loc(:), y_loc(:), x_glob(:), y_glob(:)
    complex(ep) :: ref, got
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzdotc', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 1211 + 7 * i)
        call gen_distrib_vector_z(n, mb, y_loc, y_glob, seed = 1311 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)
        call descinit_local(descy, n, 1, mb, 1, 0, 0, my_context, lld, info)

        got = (0.0_ep, 0.0_ep)
        call target_pzdotc(n, got, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1)

        if (my_rank == 0) then
            ref = zdotc(n, x_glob, 1, y_glob, 1)
            err = rel_err_scalar_z(got, ref)
            ! Complex dot has 8*n real FLOPs (4 mults + 2 adds per term).
            tol = 32.0_ep * 8.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
        end if
        deallocate(x_loc, y_loc, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzdotc
