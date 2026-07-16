program test_pdger
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dger
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gen_distrib_vector, &
                             gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdger
    implicit none

    integer, parameter :: ms(*) = [32, 80, 160]
    integer, parameter :: ns(*) = [40, 60, 120]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, locn_x, locn_y, lld_a, lld_x, lld_y
    integer :: desca(9), descx(9), descy(9)
    real(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    real(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdger', target_name, my_rank)

    alpha = 0.4_ep
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 1601 + 13 * i)
        call gen_distrib_vector(m, mb, x_loc, x_glob, seed = 1611 + 13 * i)
        call gen_distrib_vector(n, nb, y_loc, y_glob, seed = 1621 + 13 * i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locn_x = numroc_local(m, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)
        locn_y = numroc_local(n, nb, my_row, 0, my_nprow); lld_y = max(1, locn_y)

        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descx, m, 1, mb, 1, 0, 0, my_context, lld_x, info)
        call descinit_local(descy, n, 1, nb, 1, 0, 0, my_context, lld_y, info)

        call target_pdger(m, n, alpha, x_loc, 1, 1, descx, 1, &
                          y_loc, 1, 1, descy, 1, A_loc, 1, 1, desca)
        call gather_matrix(m, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n))
            A_ref = A_glob
            call dger(m, n, alpha, x_glob, 1, y_glob, 1, A_ref, m)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 32.0_ep * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdger
