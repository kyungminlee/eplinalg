program test_pdpotrf
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dpotrf
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix, &
                                scatter_matrix
    use target_scalapack, only: target_name, target_eps, target_pdpotrf
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: M_glob(:,:), dummy_loc(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    integer :: k

    call grid_init()
    call report_init('pdpotrf', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, dummy_loc, M_glob, seed = 9401 + 31*i)
        deallocate(dummy_loc)

        ! Global SPD A = M^T*M + n*I (identical on every rank since same seed).
        allocate(A_glob(n, n))
        A_glob = matmul(transpose(M_glob), M_glob)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(n, ep)
        end do

        ! Scatter A_glob into local block-cyclic slab.
        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        allocate(A_loc(max(1, locm_a), max(1, locn_a)))
        A_loc = 0.0_ep
        call scatter_matrix(n, n, mb, nb, A_glob, A_loc)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        call target_pdpotrf('U', n, A_loc, 1, 1, desca, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n))
            A_ref = A_glob
            call dpotrf('U', n, A_ref, n, info_ref)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc, A_glob, M_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpotrf
