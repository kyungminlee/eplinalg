program test_pdznrm2
    use prec_kinds,    only: ep
    use compare,       only: rel_err_scalar
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dznrm2
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z
    use target_pblas,  only: target_name, target_eps, target_pdznrm2
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9)
    complex(ep), allocatable :: x_loc(:), x_glob(:)
    real(ep) :: ref, got, err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pdznrm2', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2101 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)

        got = 0.0_ep
        call target_pdznrm2(n, got, x_loc, 1, 1, descx, 1)

        if (my_rank == 0) then
            ref = dznrm2(n, x_glob, 1)
            err = rel_err_scalar(got, ref)
            tol = 32.0_ep * 2.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
        end if
        deallocate(x_loc, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdznrm2
