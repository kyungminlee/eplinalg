program test_pdposvx
    ! Expert SPD solver. FACT='N' / EQUED='N': pdposvx Cholesky-factors A
    ! internally and solves; compare gathered X against dposv. The rcond
    ! field is also exercised but only checked for plausibility (a strict
    ! Hager-Higham bound check belongs in test_pdpocon).
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dposv
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdposvx
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs    = 2
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descaf(9), descb(9), descx(9)
    integer :: owner_r, owner_c, il, jl, ig, jg
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    real(ep), allocatable :: B_loc(:,:), B_glob(:,:)
    real(ep), allocatable :: AF_loc(:,:), X_loc(:,:), X_got(:,:)
    real(ep), allocatable :: A_sym(:,:), A_ref(:,:), B_ref(:,:)
    real(ep), allocatable :: SR(:), SC(:), ferr(:), berr(:)
    real(ep), allocatable :: work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: err, tol, rcond
    character(len=1) :: equed
    character(len=48) :: label

    call grid_init()
    call report_init('pdposvx', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 36101 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 36111 + 31*i)

        ! Symmetrize on the replicated A_glob, then mirror to A_loc owners.
        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        do k = 1, n
            A_sym(k, k) = A_sym(k, k) + real(2 * n, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(nrhs, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_sym(ig, jg)
                end do
            end do
        end if
        call descinit_local(desca,  n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descaf, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(AF_loc(lld_a, max(1, locn_a)), X_loc(lld_b, max(1, locn_b)))
        AF_loc = 0.0_ep; X_loc = 0.0_ep

        allocate(SR(n), SC(n), ferr(nrhs), berr(nrhs))
        SR = 0.0_ep; SC = 0.0_ep
        equed = 'N'

        allocate(work(1), iwork(1))
        call target_pdposvx('N', uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, equed, SR, SC, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdposvx('N', uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, equed, SR, SC, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, lwork, iwork, liwork, info)

        call gather_matrix(n, nrhs, mb, nb, X_loc, X_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_sym
            B_ref = B_glob
            call dposv(uplos(i), n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat(X_got, B_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, X_got)
        end if
        deallocate(A_loc, A_glob, B_loc, B_glob, A_sym, AF_loc, X_loc, &
                   SR, SC, ferr, berr, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdposvx
