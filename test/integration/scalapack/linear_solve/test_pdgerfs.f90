program test_pdgerfs
    ! Iterative refinement of a general linear solve. Inputs are A,
    ! the LU factor AF (from pdgetrf), B, and an initial X (from
    ! pdgetrs). pdgerfs polishes X toward the true solution using
    ! the original A and the factored AF; output X should agree
    ! with the quad reference dgesv to ~target_eps.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgesv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdgetrf, target_pdgetrs, target_pdgerfs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork
    integer :: locm_a, locn_a, locm_b, lld_a, lld_b
    integer :: desca(9), descaf(9), descb(9), descx(9)
    real(ep), allocatable :: A_loc(:,:), AF_loc(:,:), B_loc(:,:), X_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), X_got(:,:), X_ref(:,:)
    real(ep), allocatable :: ferr(:), berr(:), work(:)
    integer,  allocatable :: ipiv_got(:), ipiv_ref(:), iwork(:)
    real(ep) :: err, tol
    real(ep) :: wopt(1)
    integer  :: iwopt(1)
    character(len=48) :: label
    integer :: k, owner_r, owner_c, il, jl

    call grid_init()
    call report_init('pdgerfs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 22301 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 22311 + 31*i)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(n, ep)
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + real(n, ep)
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
        call descinit_local(desca,  n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descaf, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        ! AF holds a copy of A that we factor; A stays original for refinement.
        allocate(AF_loc(size(A_loc, 1), size(A_loc, 2)))
        AF_loc = A_loc
        allocate(ipiv_got(locm_a + mb))
        call target_pdgetrf(n, n, AF_loc, 1, 1, descaf, ipiv_got, info)

        ! Initial X = A^-1 B via pdgetrs.
        allocate(X_loc(size(B_loc, 1), size(B_loc, 2)))
        X_loc = B_loc
        call target_pdgetrs('N', n, nrhs, AF_loc, 1, 1, descaf, ipiv_got, &
                            X_loc, 1, 1, descx, info)

        allocate(ferr(nrhs), berr(nrhs))
        call target_pdgerfs('N', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, ipiv_got, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, wopt, -1, iwopt, -1, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call target_pdgerfs('N', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, ipiv_got, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, lwork, iwork, liwork, info)
        call gather_matrix(n, nrhs, mb, nb, X_loc, X_got)

        if (my_rank == 0) then
            allocate(X_ref(n, nrhs), ipiv_ref(n))
            X_ref = B_glob
            call dgesv(n, nrhs, A_glob, n, ipiv_ref, X_ref, n, info_ref)
            err = max_rel_err_mat(X_got, X_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(X_ref, ipiv_ref, X_got)
        end if
        deallocate(A_loc, AF_loc, B_loc, X_loc, A_glob, B_glob, &
                   ipiv_got, ferr, berr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgerfs
