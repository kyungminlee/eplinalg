program test_pdamax
    use prec_kinds,    only: ep
    use compare,       only: rel_err_scalar
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: idamax
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_row, numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_vector
    use target_pblas,  only: target_name, target_eps, target_pdamax
    implicit none

    integer, parameter :: cases(*) = [100, 1000, 5000]
    integer, parameter :: mb = 16
    integer :: i, n, loc_n, lld, info, indx, idx_ref
    integer :: descx(9)
    real(ep), allocatable :: x_loc(:), x_glob(:)
    real(ep) :: amax, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdamax', target_name, my_rank)

    do i = 1, size(cases)
        n = cases(i)
        call gen_distrib_vector(n, mb, x_loc, x_glob, seed = 1901 + 7 * i)

        loc_n = numroc_local(n, mb, my_row, 0, my_nprow)
        lld   = max(1, loc_n)
        call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld, info)

        amax = 0.0_ep; indx = 0
        call target_pdamax(n, amax, indx, x_loc, 1, 1, descx, 1)

        if (my_rank == 0) then
            idx_ref = idamax(n, x_glob, 1)
            ! PBLAS pXamax outputs the *signed* value at the argmax,
            ! not its absolute value (despite the docstring naming it
            ! "AMAX"). The magnitudes are what idamax compared to find
            ! the index; the value carried back is x_glob(indx) itself.
            !
            ! Tie-tolerant index check: distinct indices are accepted
            ! when the magnitudes coincide (continuous RNG makes this
            ! astronomically unlikely, but it's the right rule).
            if (indx == idx_ref) then
                err = rel_err_scalar(amax, x_glob(indx))
            else if (abs(abs(x_glob(indx)) - abs(x_glob(idx_ref))) <= &
                     2.0_ep * target_eps * abs(x_glob(idx_ref))) then
                err = rel_err_scalar(amax, x_glob(indx))
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
end program test_pdamax
