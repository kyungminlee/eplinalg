program test_pzgeru
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zgeru
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzgeru
    implicit none

    integer, parameter :: ms(*) = [24, 64, 100]
    integer, parameter :: ns(*) = [32, 48, 120]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, locn_x, locn_y, lld_a, lld_x, lld_y
    integer :: desca(9), descx(9), descy(9)
    complex(ep), allocatable :: A_loc(:,:), x_loc(:), y_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), y_glob(:), &
                                A_got(:,:), A_ref(:,:)
    complex(ep) :: alpha
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzgeru', target_name, my_rank)

    alpha = cmplx(0.4_ep, 0.2_ep, ep)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix_z(m, n, mb, nb, A_loc, A_glob, seed = 5601 + 23 * i)
        call gen_distrib_vector_z(m, mb, x_loc, x_glob, seed = 5611 + 23 * i)
        call gen_distrib_vector_z(n, nb, y_loc, y_glob, seed = 5621 + 23 * i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locn_x = numroc_local(m, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)
        locn_y = numroc_local(n, nb, my_row, 0, my_nprow); lld_y = max(1, locn_y)

        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descx, m, 1, mb, 1, 0, 0, my_context, lld_x, info)
        call descinit_local(descy, n, 1, nb, 1, 0, 0, my_context, lld_y, info)

        call target_pzgeru(m, n, alpha, x_loc, 1, 1, descx, 1, &
                           y_loc, 1, 1, descy, 1, A_loc, 1, 1, desca)
        call gather_matrix_z(m, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n))
            A_ref = A_glob
            call zgeru(m, n, alpha, x_glob, 1, y_glob, 1, A_ref, m)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 32.0_ep * 8.0_ep * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc, x_loc, y_loc, A_glob, x_glob, y_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgeru
