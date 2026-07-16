program test_pzposv
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zposv
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzposv
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs = 3
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, k
    integer :: locm_a, locn_a, lld_a, locm_b, lld_b
    integer :: desca(9), descb(9)
    integer :: owner_r, owner_c, il, jl, ig, jg
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), B_got(:,:)
    complex(ep), allocatable :: A_herm(:,:), A_ref(:,:), B_ref(:,:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pzposv', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n,    mb, nb, A_loc, A_glob, seed = 12301 + 31*i)
        call gen_distrib_matrix_z(n, nrhs, mb, nb, B_loc, B_glob, seed = 12311 + 31*i)

        ! Hermitize and add diagonal boost for HPD.
        allocate(A_herm(n, n))
        A_herm = 0.5_ep * (A_glob + conjg(transpose(A_glob)))
        do k = 1, n
            A_herm(k, k) = A_herm(k, k) + cmplx(real(2 * n, ep), 0.0_ep, kind=ep)
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
                    if (owner_r == my_row) A_loc(il, jl) = A_herm(ig, jg)
                end do
            end do
        end if
        call descinit_local(desca, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        call target_pzposv(uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                           B_loc, 1, 1, descb, info)
        call gather_matrix_z(n, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_herm
            B_ref = B_glob
            call zposv(uplos(i), n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,a,a,i0,a,i0)') 'uplo=', uplos(i), &
                    ',n=', n, ',nrhs=', nrhs
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, B_ref, B_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob, A_herm)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzposv
