program test_pzhegst
    ! Hermitian generalized eigenvalue reduction (complex mirror of
    ! pdsygst): given Hermitian A and Cholesky factor of HPD B (in B),
    ! reduce A to standard form. Compare gathered A against zhegst.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: zpotrf, zhegst
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pzpotrf, target_pzhegst
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: ibtypes(*) = [1, 2, 3]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), descb(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), A_got(:,:)
    complex(ep), allocatable :: A_herm(:,:), B_herm(:,:), A_ref(:,:), B_ref(:,:)
    real(ep) :: err, tol, scale
    character(len=48) :: label

    call grid_init()
    call report_init('pzhegst', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 25101 + 31*i)
        call gen_distrib_matrix_z(n, n, mb, nb, B_loc, B_glob, seed = 25111 + 31*i)

        allocate(A_herm(n, n), B_herm(n, n))
        A_herm = 0.5_ep * (A_glob + conjg(transpose(A_glob)))
        B_herm = 0.5_ep * (B_glob + conjg(transpose(B_glob)))
        do k = 1, n
            B_herm(k, k) = B_herm(k, k) + cmplx(real(2 * n, ep), 0.0_ep, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) then
                        A_loc(il, jl) = A_herm(ig, jg)
                        B_loc(il, jl) = B_herm(ig, jg)
                    end if
                end do
            end do
        end if
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        call target_pzpotrf(uplos(i), n, B_loc, 1, 1, descb, info)
        call target_pzhegst(ibtypes(i), uplos(i), n, A_loc, 1, 1, desca, &
                            B_loc, 1, 1, descb, scale, info)
        call gather_matrix_z(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, n))
            A_ref = A_herm
            B_ref = B_herm
            call zpotrf(uplos(i), n, B_ref, n, info_ref)
            call zhegst(ibtypes(i), uplos(i), n, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,a,a,i0)') 'ibtype=', ibtypes(i), &
                ',uplo=', uplos(i), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, A_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob, A_herm, B_herm)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzhegst
