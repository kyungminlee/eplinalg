program test_pdlacpy
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdlacpy
    implicit none

    integer, parameter :: ms(*) = [32, 64, 96]
    integer, parameter :: ns(*) = [40, 48, 80]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:), A_glob(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdlacpy', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 11101 + 31*i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, m, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(B_loc(max(1, locm_a), max(1, locn_a)))
        B_loc = 0.0_ep

        call target_pdlacpy('A', m, n, A_loc, 1, 1, desca, B_loc, 1, 1, descb)
        call gather_matrix(m, n, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            err = max_rel_err_mat(B_got, A_glob)
            tol = target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(B_got)
        end if
        deallocate(A_loc, B_loc, A_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdlacpy
