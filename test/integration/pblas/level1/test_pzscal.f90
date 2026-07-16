program test_pzscal
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zscal
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z, gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzscal
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info
    integer :: descx(9)
    complex(ep), allocatable :: x_loc(:), x_glob(:), x_got(:), x_ref(:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=32) :: label

    call grid_init()
    call report_init('pzscal', target_name, my_rank)

    alpha = cmplx(0.75_ep, -0.25_ep, ep)
    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2501 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)

        call target_pzscal(n, alpha, x_loc, 1, 1, descx, 1)
        call gather_vector_z(n, mb, x_loc, x_got)

        if (my_rank == 0) then
            allocate(x_ref(n))
            x_ref = x_glob
            call zscal(n, alpha, x_ref, 1)
            err = max_rel_err_vec_z(x_got, x_ref)
            ! One complex multiply per element: 4 mults + 2 adds = 8 FLOPs.
            tol = 32.0_ep * 8.0_ep * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(x_ref, x_got)
        end if
        deallocate(x_loc, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzscal
