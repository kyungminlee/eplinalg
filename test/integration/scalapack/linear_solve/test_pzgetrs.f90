program test_pzgetrs
    ! Z-mirror of test_pdgetrs.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgetrf, zgetrs
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, &
                                target_pzgetrf, target_pzgetrs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs  = 2
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref
    integer :: locm_a, locm_b, lld_a, lld_b
    integer :: desca(9), descb(9)
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), B_got(:,:), B_ref(:,:)
    integer,     allocatable :: ipiv_got(:), ipiv_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label
    integer :: k, owner_r, owner_c, il, jl

    call grid_init()
    call report_init('pzgetrs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n,    mb, nb, A_loc, A_glob, seed = 21801 + 31*i)
        call gen_distrib_matrix_z(n, nrhs, mb, nb, B_loc, B_glob, seed = 21811 + 31*i)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + cmplx(real(n, ep), 0.0_ep, ep)
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + cmplx(real(n, ep), 0.0_ep, ep)
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
        call descinit_local(desca, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(ipiv_got(locm_a + mb))
        call target_pzgetrf(n, n, A_loc, 1, 1, desca, ipiv_got, info)
        call target_pzgetrs('N', n, nrhs, A_loc, 1, 1, desca, ipiv_got, &
                            B_loc, 1, 1, descb, info)
        call gather_matrix_z(n, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(B_ref(n, nrhs), ipiv_ref(n))
            B_ref = B_glob
            call zgetrf(n, n, A_glob, n, ipiv_ref, info_ref)
            call zgetrs('N', n, nrhs, A_glob, n, ipiv_ref, B_ref, n, info_ref)
            err = max_rel_err_mat_z(B_got, B_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(B_ref, ipiv_ref, B_got)
        end if
        deallocate(A_loc, B_loc, A_glob, B_glob, ipiv_got)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgetrs
