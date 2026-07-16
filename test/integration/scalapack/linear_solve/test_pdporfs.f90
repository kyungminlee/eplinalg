program test_pdporfs
    ! Iterative refinement of a Cholesky-based SPD solve.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dposv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdpotrf, target_pdpotrs, target_pdporfs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork
    integer :: locm_a, locn_a, locm_b, lld_a, lld_b
    integer :: desca(9), descaf(9), descb(9), descx(9)
    real(ep), allocatable :: A_loc(:,:), AF_loc(:,:), B_loc(:,:), X_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), M_glob(:,:), dummy(:,:)
    real(ep), allocatable :: X_got(:,:), X_ref(:,:)
    real(ep), allocatable :: ferr(:), berr(:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: err, tol
    real(ep) :: wopt(1)
    integer  :: iwopt(1)
    character(len=48) :: label
    integer :: ig, jg, owner_r, owner_c, il, jl, k

    call grid_init()
    call report_init('pdporfs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, dummy, M_glob, seed = 22501 + 31*i)
        deallocate(dummy)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 22511 + 31*i)
        allocate(A_glob(n, n))
        A_glob = matmul(transpose(M_glob), M_glob)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(n, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
        allocate(A_loc(max(1, locm_a), max(1, locn_a)))
        A_loc = 0.0_ep
        do jg = 1, n
            call g2l(jg, nb, my_npcol, owner_c, jl)
            if (owner_c /= my_col) cycle
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                if (owner_r == my_row) A_loc(il, jl) = A_glob(ig, jg)
            end do
        end do
        call descinit_local(desca,  n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descaf, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(AF_loc(size(A_loc, 1), size(A_loc, 2)))
        AF_loc = A_loc
        call target_pdpotrf('U', n, AF_loc, 1, 1, descaf, info)

        allocate(X_loc(size(B_loc, 1), size(B_loc, 2)))
        X_loc = B_loc
        call target_pdpotrs('U', n, nrhs, AF_loc, 1, 1, descaf, &
                            X_loc, 1, 1, descx, info)

        allocate(ferr(nrhs), berr(nrhs))
        call target_pdporfs('U', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, wopt, -1, iwopt, -1, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call target_pdporfs('U', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, lwork, iwork, liwork, info)
        call gather_matrix(n, nrhs, mb, nb, X_loc, X_got)

        if (my_rank == 0) then
            allocate(X_ref(n, nrhs))
            X_ref = B_glob
            call dposv('U', n, nrhs, A_glob, n, X_ref, n, info_ref)
            err = max_rel_err_mat(X_got, X_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(X_ref, X_got)
        end if
        deallocate(A_loc, AF_loc, B_loc, X_loc, A_glob, B_glob, M_glob, &
                   ferr, berr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdporfs
