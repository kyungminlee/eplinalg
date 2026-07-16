program test_pdpocon
    ! 1-norm reciprocal condition-number estimator for SPD matrices.
    ! pdpocon needs a Cholesky factor in A; we factor first via the
    ! existing target_pdpotrf wrapper, then call target_pdpocon. Same
    ! Hager-Higham bound semantics as test_pdgecon.
    use prec_kinds,        only: ep
    use compare,           only: rel_err_scalar
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use cond_helpers,      only: true_kappa1_posdef
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix
    use target_scalapack,  only: target_name, target_pdpocon, target_pdpotrf, &
                                 target_pdlansy
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_sym(:,:), A_for_kappa(:,:)
    real(ep), allocatable :: work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: anorm, rcond, kappa_est, kappa_true, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdpocon', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 47201 + 31*i)
        ! Symmetrize A_glob (replicated on every rank), boost diagonal.
        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        do k = 1, n
            A_sym(k, k) = A_sym(k, k) + real(2 * n, ep)
        end do
        ! Reflect into A_loc owners.
        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
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
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(work(max(1, locn_a)))
        anorm = target_pdlansy('1', uplos(i), n, A_loc, 1, 1, desca, work)
        deallocate(work)

        ! Cholesky-factor in place — pdpocon expects the factor.
        call target_pdpotrf(uplos(i), n, A_loc, 1, 1, desca, info)
        if (info /= 0 .and. my_rank == 0) then
            write(*, '(a,i0,a,i0)') 'pdpotrf failed n=', n, ' info=', info
        end if

        allocate(work(1), iwork(1))
        call target_pdpocon(uplos(i), n, A_loc, 1, 1, desca, anorm, rcond, &
                            work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdpocon(uplos(i), n, A_loc, 1, 1, desca, anorm, rcond, &
                            work, lwork, iwork, liwork, info)
        kappa_est = 1.0_ep / rcond

        if (my_rank == 0) then
            allocate(A_for_kappa(n, n))
            A_for_kappa = A_sym
            call true_kappa1_posdef(uplos(i), n, A_for_kappa, n, kappa_true, info_ref)
            err = rel_err_scalar(kappa_est, kappa_true)
            tol = 0.7_ep
            write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_for_kappa)
        end if
        deallocate(A_loc, A_glob, A_sym, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpocon
