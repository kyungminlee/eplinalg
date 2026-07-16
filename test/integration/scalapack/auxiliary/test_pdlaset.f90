program test_pdlaset
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdlaset
    implicit none

    integer, parameter :: ms(*) = [32, 64, 96]
    integer, parameter :: ns(*) = [40, 48, 80]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_got(:,:), A_ref(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label
    integer :: k

    call grid_init()
    call report_init('pdlaset', target_name, my_rank)

    alpha = 0.125_ep; beta = 7.0_ep

    do i = 1, size(ms)
        m = ms(i); n = ns(i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(A_loc(max(1, locm_a), max(1, locn_a)))
        A_loc = -99.0_ep

        call target_pdlaset('A', m, n, alpha, beta, A_loc, 1, 1, desca)
        call gather_matrix(m, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n))
            A_ref = alpha
            do k = 1, min(m, n)
                A_ref(k, k) = beta
            end do
            err = max_rel_err_mat(A_got, A_ref)
            tol = target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdlaset
