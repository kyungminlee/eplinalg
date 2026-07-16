program test_pdgesvx
    ! Expert general solver. FACT='N' / EQUED='N' / TRANS='N': pdgesvx
    ! factors A internally and solves; compare gathered X against dgesv.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dgesv
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdgesvx
    implicit none

    integer, parameter :: ns(*)   = [32, 64, 96]
    integer, parameter :: nrhs    = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descaf(9), descb(9), descx(9)
    integer :: owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    real(ep), allocatable :: B_loc(:,:), B_glob(:,:)
    real(ep), allocatable :: AF_loc(:,:), X_loc(:,:), X_got(:,:)
    real(ep), allocatable :: A_ref(:,:), B_ref(:,:)
    real(ep), allocatable :: R_vec(:), C_vec(:), ferr(:), berr(:)
    real(ep), allocatable :: work(:)
    integer,  allocatable :: ipiv(:), iwork(:), ipiv_ref(:)
    real(ep) :: err, tol, rcond
    character(len=1) :: equed
    character(len=48) :: label

    call grid_init()
    call report_init('pdgesvx', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 29101 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 29111 + 31*i)

        ! Diagonally-bias A on every rank (A_glob is replicated by gen_distrib_matrix)
        ! and on the owner rank's local copy.
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(2 * n, ep)
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + real(2 * n, ep)
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(nrhs, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        call descinit_local(desca,  n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descaf, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(AF_loc(lld_a, max(1, locn_a)), X_loc(lld_b, max(1, locn_b)))
        AF_loc = 0.0_ep; X_loc = 0.0_ep

        allocate(ipiv(lld_a + nb), R_vec(n), C_vec(n), &
                 ferr(nrhs), berr(nrhs))
        R_vec = 0.0_ep; C_vec = 0.0_ep
        equed = 'N'

        allocate(work(1), iwork(1))
        call target_pdgesvx('N', 'N', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, ipiv, equed, R_vec, C_vec, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdgesvx('N', 'N', n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, ipiv, equed, R_vec, C_vec, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, lwork, iwork, liwork, info)

        call gather_matrix(n, nrhs, mb, nb, X_loc, X_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs), ipiv_ref(n))
            A_ref = A_glob; B_ref = B_glob
            call dgesv(n, nrhs, A_ref, n, ipiv_ref, B_ref, n, info_ref)
            err = max_rel_err_mat(X_got, B_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, ipiv_ref, X_got)
        end if
        deallocate(A_loc, A_glob, B_loc, B_glob, AF_loc, X_loc, &
                   ipiv, R_vec, C_vec, ferr, berr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgesvx
