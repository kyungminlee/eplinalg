program test_pzamax
    use prec_kinds,    only: ep
    use compare,       only: rel_err_scalar_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: izamax
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector_z
    use target_pblas,  only: target_name, target_eps, target_pzamax
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info, indx, idx_ref
    integer :: descx(9)
    complex(ep), allocatable :: x_loc(:), x_glob(:)
    complex(ep) :: amax
    real(ep) :: mag_indx, mag_ref, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzamax', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector_z(n, mb, x_loc, x_glob, seed = 2701 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)

        amax = (0.0_ep, 0.0_ep); indx = 0
        call target_pzamax(n, amax, indx, x_loc, 1, 1, descx, 1)

        if (my_rank == 0) then
            idx_ref = izamax(n, x_glob, 1)
            ! BLAS convention: |x| = |re| + |im| (1-norm), not Euclidean.
            mag_indx = abs(real(x_glob(indx), ep)) + abs(aimag(x_glob(indx)))
            mag_ref  = abs(real(x_glob(idx_ref), ep)) + abs(aimag(x_glob(idx_ref)))
            if (indx == idx_ref) then
                err = rel_err_scalar_z(amax, x_glob(indx))
            else if (abs(mag_indx - mag_ref) <= 2.0_ep * target_eps * mag_ref) then
                err = rel_err_scalar_z(amax, x_glob(indx))
            else
                err = huge(1.0_ep)
            end if
            tol = 2.0_ep * target_eps
            write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',indx=', indx, &
                ',ref=', idx_ref
            call report_case(trim(label), err, tol)
        end if
        deallocate(x_loc, x_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzamax
