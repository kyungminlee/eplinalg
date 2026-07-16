program test_pdtrtrs
    ! Triangular solve A*X = B (UPLO='U', DIAG='N', TRANS='N').
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dtrtrs
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdtrtrs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref
    integer :: locm_a, locn_a, locm_b, lld_a, lld_b
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), B_got(:,:), B_ref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label
    integer :: k, owner_r, owner_c, il, jl, ig, jg

    call grid_init()
    call report_init('pdtrtrs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 23001 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 23011 + 31*i)
        ! Make A upper-triangular and well-conditioned: zero strict lower,
        ! boost diagonal.
        do jg = 1, n
            do ig = jg + 1, n
                A_glob(ig, jg) = 0.0_ep
                call g2l(ig, mb, my_nprow, owner_r, il)
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_r == my_row .and. owner_c == my_col) A_loc(il, jl) = 0.0_ep
            end do
        end do
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(n, ep)
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + real(n, ep)
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
        call descinit_local(desca, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        call target_pdtrtrs('U', 'N', 'N', n, nrhs, A_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, info)
        call gather_matrix(n, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(B_ref(n, nrhs))
            B_ref = B_glob
            call dtrtrs('U', 'N', 'N', n, nrhs, A_glob, n, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_ref, B_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtrtrs
