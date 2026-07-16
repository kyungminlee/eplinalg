program test_pdgehd2
    ! Hessenberg reduction. Compare gathered post-call A (encodes H +
    ! reflectors) against dgehd2.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgehrd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgehd2
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, ilo, ihi, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: tau_got(:), tau_ref(:), work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdgehd2', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i); ilo = 1; ihi = n
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 15201 + 31*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        ! pdgehd2's documented LWORK >= NB + max(NpA0, NB) under-states
        ! actual usage on small grids. Allocate a generous buffer.
        lwork = max(1024, n * nb)
        allocate(tau_got(max(1, locn_a)), work(lwork))
        call target_pdgehd2(n, ilo, ihi, A_loc, 1, 1, desca, tau_got, &
                            work, lwork, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), tau_ref(max(1, n - 1)), &
                     work_ref(max(1, n * 64)))
            A_ref = A_glob
            call dgehrd(n, ilo, ihi, A_ref, n, tau_ref, work_ref, &
                        size(work_ref), info_ref)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, tau_ref, work_ref, A_got)
        end if
        deallocate(A_loc, A_glob, tau_got, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgehd2
