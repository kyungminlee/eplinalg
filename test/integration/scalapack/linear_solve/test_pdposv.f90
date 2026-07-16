program test_pdposv
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dposv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdposv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 3
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, k
    integer :: locm_a, locn_a, lld_a, locm_b, lld_b
    integer :: desca(9), descb(9)
    integer :: owner_r, owner_c, il, jl, ig, jg
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), B_got(:,:)
    real(ep), allocatable :: A_sym(:,:), A_ref(:,:), B_ref(:,:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pdposv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n,    mb, nb, A_loc, A_glob, seed = 12001 + 31*i)
        call gen_distrib_matrix(n, nrhs, mb, nb, B_loc, B_glob, seed = 12011 + 31*i)

        ! Symmetrize and add a diagonal boost so A is SPD.
        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))
        do k = 1, n
            A_sym(k, k) = A_sym(k, k) + real(2 * n, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
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
        call descinit_local(desca, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        call target_pdposv(uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                           B_loc, 1, 1, descb, info)
        call gather_matrix(n, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_sym
            B_ref = B_glob
            call dposv(uplos(i), n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
                    ',n=', n, ',nrhs=', nrhs
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, B_ref, B_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob, A_sym)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdposv
