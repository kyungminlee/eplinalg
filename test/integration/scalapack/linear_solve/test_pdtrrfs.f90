program test_pdtrrfs
    ! Iterative-refinement error bounds for triangular system.
    ! Construct a triangular T and exact X, set B = T*X, then call
    ! pdtrrfs to compute FERR/BERR. Compare against dtrrfs by running
    ! both with the same (T, X, B) and checking that the two BERR
    ! vectors agree to within a generous factor (the Hager estimator
    ! has algorithmic slop, so this is a smoke + sanity test rather
    ! than a tight precision check).
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dtrrfs
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdtrrfs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, k, j
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descb(9), descx(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: T_loc(:,:), T_glob(:,:)
    real(ep), allocatable :: X_loc(:,:), X_glob(:,:)
    real(ep), allocatable :: B_loc(:,:), B_glob(:,:)
    real(ep), allocatable :: ferr(:), berr(:), ferr_ref(:), berr_ref(:)
    real(ep), allocatable :: T_ref(:,:), X_ref(:,:), B_ref(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: iwork(:), iwork_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdtrrfs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, T_loc, T_glob, seed = 30101 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, X_loc, X_glob, seed = 30111 + 31*i)

        ! Force upper-triangular T globally and locally; bias diagonal.
        do jg = 1, n
            do ig = jg + 1, n
                T_glob(ig, jg) = 0.0_ep
            end do
            T_glob(jg, jg) = T_glob(jg, jg) + real(2 * n, ep)
        end do
        if (size(T_loc, 1) > 0 .and. size(T_loc, 2) > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) T_loc(il, jl) = T_glob(ig, jg)
                end do
            end do
        end if

        ! Compute B = T*X globally so the system is exact, then distribute.
        allocate(B_glob(n, nrhs))
        B_glob = matmul(T_glob, X_glob)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(nrhs, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        allocate(B_loc(lld_b, max(1, locn_b)))
        B_loc = 0.0_ep
        do jg = 1, nrhs
            call g2l(jg, nb, my_npcol, owner_c, jl)
            if (owner_c /= my_col) cycle
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                if (owner_r == my_row) B_loc(il, jl) = B_glob(ig, jg)
            end do
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(ferr(nrhs), berr(nrhs))
        allocate(work(1), iwork(1))
        call target_pdtrrfs('U', 'N', 'N', n, nrhs, T_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdtrrfs('U', 'N', 'N', n, nrhs, T_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, lwork, iwork, liwork, info)

        if (my_rank == 0) then
            allocate(T_ref(n, n), X_ref(n, nrhs), B_ref(n, nrhs), &
                     ferr_ref(nrhs), berr_ref(nrhs), &
                     work_ref(3 * n), iwork_ref(n))
            T_ref = T_glob; X_ref = X_glob; B_ref = B_glob
            call dtrrfs('U', 'N', 'N', n, nrhs, T_ref, n, &
                        B_ref, n, X_ref, n, ferr_ref, berr_ref, &
                        work_ref, iwork_ref, info_ref)
            ! BERR is what's most stable here. Compare absolute logarithm:
            ! both should be ~unit roundoff since X is exact. We accept a
            ! constant factor (10^3) — Hager estimator + parallel reduction
            ! reorder gives appreciable variation.
            err = 0.0_ep
            do j = 1, nrhs
                if (berr(j) /= 0.0_ep .and. berr_ref(j) /= 0.0_ep) then
                    err = max(err, abs(log(berr(j) / berr_ref(j))) / log(10.0_ep))
                end if
            end do
            ! err is now |orders of magnitude difference|.
            tol = 3.0_ep  ! allow up to 3 orders of magnitude (Hager slop)
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(T_ref, X_ref, B_ref, ferr_ref, berr_ref, work_ref, iwork_ref)
        end if
        deallocate(T_loc, T_glob, X_loc, X_glob, B_loc, B_glob, &
                   ferr, berr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtrrfs
