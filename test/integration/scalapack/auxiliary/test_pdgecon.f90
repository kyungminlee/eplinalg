program test_pdgecon
    ! 1-norm reciprocal condition-number estimator. PDGECON returns a
    ! Hager-Higham underestimate: 1/κ_true ≤ rcond_est ≤ 3/κ_true (i.e.
    ! κ_true/3 ≤ κ_est ≤ κ_true). κ_true is computed exactly at quad
    ! precision via reference dlange + dgetrf + dgetri on rank 0.
    use prec_kinds,        only: ep
    use compare,           only: rel_err_scalar
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use cond_helpers,      only: true_kappa1_general
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix
    use target_scalapack,  only: target_name, target_pdgecon, target_pdgetrf, &
                                 target_pdlange
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_for_kappa(:,:)
    real(ep), allocatable :: work(:)
    integer,  allocatable :: ipiv(:), iwork(:)
    real(ep) :: anorm, rcond, kappa_est, kappa_true, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdgecon', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        ! Replicated A_glob on every rank; A_loc is the distribution.
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 47101 + 31*i)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(2 * n, ep)
        end do
        ! Mirror diagonal boost into A_loc owners.
        block
            integer :: ig, owner_r, owner_c, il, jl
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                call g2l(ig, nb, my_npcol, owner_c, jl)
                if (owner_r == my_row .and. owner_c == my_col) then
                    A_loc(il, jl) = A_loc(il, jl) + real(2 * n, ep)
                end if
            end do
        end block

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        ! ANORM via the existing pdlange wrapper, taken on the
        ! original A *before* the in-place LU.
        allocate(work(max(1, locn_a)))
        anorm = target_pdlange('1', n, n, A_loc, 1, 1, desca, work)
        deallocate(work)

        ! pdgecon expects A to hold the LU factor (output of pdgetrf).
        allocate(ipiv(lld_a + nb))
        call target_pdgetrf(n, n, A_loc, 1, 1, desca, ipiv, info)
        if (info /= 0 .and. my_rank == 0) then
            write(*, '(a,i0,a,i0)') 'pdgetrf failed n=', n, ' info=', info
        end if

        ! Workspace query, then call.
        allocate(work(1), iwork(1))
        call target_pdgecon('1', n, A_loc, 1, 1, desca, anorm, rcond, &
                            work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdgecon('1', n, A_loc, 1, 1, desca, anorm, rcond, &
                            work, lwork, iwork, liwork, info)
        kappa_est = 1.0_ep / rcond

        if (my_rank == 0) then
            allocate(A_for_kappa(n, n))
            A_for_kappa = A_glob
            call true_kappa1_general(n, A_for_kappa, n, kappa_true, info_ref)
            ! LACON underestimates κ by up to a factor of 3:
            !   κ_true / 3 ≤ κ_est ≤ κ_true
            ! Report rel-err vs κ_true; tolerance = 0.7 (covers the
            ! worst-case 0.333 ratio with margin for fp noise).
            err = rel_err_scalar(kappa_est, kappa_true)
            tol = 0.7_ep
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_for_kappa)
        end if
        deallocate(A_loc, A_glob, work, iwork, ipiv)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgecon
