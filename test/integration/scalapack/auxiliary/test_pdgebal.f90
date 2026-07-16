program test_pdgebal
    ! Eigenvalue-prep balancing of a general matrix. Output: ILO, IHI,
    ! SCALE vector, and a balanced A. Compare against serial dgebal.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat, max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgebal
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgebal
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, ilo, ihi, ilo_ref, ihi_ref
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: scale_got(:), scale_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdgebal', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 22901 + 31*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(scale_got(n))
        call target_pdgebal('B', n, A_loc, desca, ilo, ihi, scale_got, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), scale_ref(n))
            A_ref = A_glob
            call dgebal('B', n, A_ref, n, ilo_ref, ihi_ref, scale_ref, info_ref)
            err = max(max_rel_err_mat(A_got, A_ref), &
                      max_rel_err_vec(scale_got, scale_ref))
            if (ilo /= ilo_ref) err = max(err, 1.0_ep)
            if (ihi /= ihi_ref) err = max(err, 1.0_ep)
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, scale_ref, A_got)
        end if
        deallocate(A_loc, A_glob, scale_got)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgebal
