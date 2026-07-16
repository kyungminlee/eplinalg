program test_pztranu
    ! pztranu: C := beta*C + alpha*A^T (un-conjugated transpose).
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pztranu
    implicit none

    integer, parameter :: ms(*) = [32, 64, 96]
    integer, parameter :: ns(*) = [40, 48, 80]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, locm_c, locn_c, lld_a, lld_c
    integer :: desca(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pztranu', target_name, my_rank)

    alpha = cmplx(0.7_ep,  0.2_ep, kind=ep)
    beta  = cmplx(0.3_ep, -0.1_ep, kind=ep)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix_z(n, m, mb, nb, A_loc, A_glob, seed = 17501 + 29*i)
        call gen_distrib_matrix_z(m, n, mb, nb, C_loc, C0,    seed = 17511 + 29*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(m, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_c = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(n, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
        call descinit_local(desca, n, m, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descc, m, n, mb, nb, 0, 0, my_context, lld_c, info)

        call target_pztranu(m, n, alpha, A_loc, 1, 1, desca, &
                            beta, C_loc, 1, 1, descc)
        call gather_matrix_z(m, n, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(m, n))
            C_ref = beta * C0 + alpha * transpose(A_glob)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 16.0_ep * real(max(m, n), ep) * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, C_loc, A_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pztranu
