program test_pztrrfs
    ! Iterative-refinement error bounds for complex triangular system.
    ! See test_pdtrrfs for design.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: ztrrfs
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z
    use target_scalapack,  only: target_name, target_eps, target_pztrrfs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, lrwork, j
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descb(9), descx(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: T_loc(:,:), T_glob(:,:)
    complex(ep), allocatable :: X_loc(:,:), X_glob(:,:)
    complex(ep), allocatable :: B_loc(:,:), B_glob(:,:)
    real(ep), allocatable :: ferr(:), berr(:), ferr_ref(:), berr_ref(:)
    complex(ep), allocatable :: T_ref(:,:), X_ref(:,:), B_ref(:,:)
    complex(ep), allocatable :: work(:), work_ref(:)
    real(ep), allocatable :: rwork(:), rwork_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pztrrfs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n,    mb, nb, T_loc, T_glob, seed = 30201 + 31*i)
        call gen_distrib_matrix_z(n, nrhs, mb, nb, X_loc, X_glob, seed = 30211 + 31*i)

        do jg = 1, n
            do ig = jg + 1, n
                T_glob(ig, jg) = (0.0_ep, 0.0_ep)
            end do
            T_glob(jg, jg) = T_glob(jg, jg) + cmplx(real(2 * n, ep), 0.0_ep, ep)
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

        allocate(B_glob(n, nrhs))
        B_glob = matmul(T_glob, X_glob)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(nrhs, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        allocate(B_loc(lld_b, max(1, locn_b)))
        B_loc = (0.0_ep, 0.0_ep)
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
        allocate(work(1), rwork(1))
        call target_pztrrfs('U', 'N', 'N', n, nrhs, T_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, -1, rwork, -1, info)
        lwork  = max(1, int(real(work(1))))
        lrwork = max(1, int(rwork(1)))
        deallocate(work, rwork)
        allocate(work(lwork), rwork(lrwork))
        call target_pztrrfs('U', 'N', 'N', n, nrhs, T_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            ferr, berr, work, lwork, rwork, lrwork, info)

        if (my_rank == 0) then
            allocate(T_ref(n, n), X_ref(n, nrhs), B_ref(n, nrhs), &
                     ferr_ref(nrhs), berr_ref(nrhs), &
                     work_ref(2 * n), rwork_ref(n))
            T_ref = T_glob; X_ref = X_glob; B_ref = B_glob
            call ztrrfs('U', 'N', 'N', n, nrhs, T_ref, n, &
                        B_ref, n, X_ref, n, ferr_ref, berr_ref, &
                        work_ref, rwork_ref, info_ref)
            err = 0.0_ep
            do j = 1, nrhs
                if (berr(j) /= 0.0_ep .and. berr_ref(j) /= 0.0_ep) then
                    err = max(err, abs(log(berr(j) / berr_ref(j))) / log(10.0_ep))
                end if
            end do
            tol = 3.0_ep
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(T_ref, X_ref, B_ref, ferr_ref, berr_ref, work_ref, rwork_ref)
        end if
        deallocate(T_loc, T_glob, X_loc, X_glob, B_loc, B_glob, &
                   ferr, berr, work, rwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrrfs
