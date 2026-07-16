program test_pzdscal
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zdscal
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z, gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzdscal
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9)
    complex(ep), allocatable :: x_loc(:), x_glob(:), x_got(:), x_ref(:)
    real(ep) :: alpha, err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzdscal', target_name, my_rank)

    alpha = 0.75_ep
    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2601 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)

        call target_pzdscal(n, alpha, x_loc, 1, 1, descx, 1)
        call gather_vector_z(n, mb, x_loc, x_got)

        if (my_rank == 0) then
            allocate(x_ref(n))
            x_ref = x_glob
            call zdscal(n, alpha, x_ref, 1)
            err = max_rel_err_vec_z(x_got, x_ref)
            ! Real-by-complex scale: 2 mults per element.
            tol = 32.0_ep * 2.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, x_got)
        end if
        deallocate(x_loc, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzdscal
