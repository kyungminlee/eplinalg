program test_pdpotri
    ! Cholesky factor + invert. Compare A^{-1} against dpotrf+dpotri.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dpotrf, dpotri
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdpotrf, target_pdpotri
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: owner_r, owner_c, il, jl, ig, jg
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:)
    real(ep), allocatable :: A_sym(:,:), A_ref(:,:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pdpotri', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 12101 + 31*i)

        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        do k = 1, n
            A_sym(k, k) = A_sym(k, k) + real(2 * n, ep)
        end do

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

        call target_pdpotrf(uplos(i), n, A_loc, 1, 1, desca, info)
        call target_pdpotri(uplos(i), n, A_loc, 1, 1, desca, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n))
            A_ref = A_sym
            call dpotrf(uplos(i), n, A_ref, n, info_ref)
            call dpotri(uplos(i), n, A_ref, n, info_ref)
            ! Only the named triangle is meaningful; mask off the other.
            if (uplos(i) == 'U') then
                do jg = 1, n
                    do ig = jg + 1, n
                        A_ref(ig, jg) = 0.0_ep
                        A_got(ig, jg) = 0.0_ep
                    end do
                end do
            else
                do jg = 1, n
                    do ig = 1, jg - 1
                        A_ref(ig, jg) = 0.0_ep
                        A_got(ig, jg) = 0.0_ep
                    end do
                end do
            end if
            err = max_rel_err_mat(A_got, A_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, A_got)
        end if
        deallocate(A_loc, A_glob, A_sym)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdpotri
