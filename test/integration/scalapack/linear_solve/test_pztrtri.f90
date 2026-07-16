program test_pztrtri
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: ztrtri
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pztrtri
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = [character(len=1) :: 'U', 'L']
    character(len=1), parameter :: diags(*) = [character(len=1) :: 'N', 'U']
    integer :: i, ku, kd, n, info, info_ref, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    complex(ep) :: boost
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pztrtri', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        do ku = 1, size(uplos)
            do kd = 1, size(diags)
                call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, &
                                          seed = 16201 + 13*i + 17*ku + 19*kd)
                boost = cmplx(real(2 * n, ep), 0.0_ep, kind=ep)
                do k = 1, n
                    A_glob(k, k) = A_glob(k, k) + boost
                    call g2l(k, mb, my_nprow, owner_r, il)
                    call g2l(k, nb, my_npcol, owner_c, jl)
                    if (owner_r == my_row .and. owner_c == my_col) then
                        A_loc(il, jl) = A_loc(il, jl) + boost
                    end if
                end do

                locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
                locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
                call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

                call target_pztrtri(uplos(ku), diags(kd), n, A_loc, 1, 1, desca, info)
                call gather_matrix_z(n, n, mb, nb, A_loc, A_got)

                if (my_rank == 0) then
                    allocate(A_ref(n, n))
                    A_ref = A_glob
                    call ztrtri(uplos(ku), diags(kd), n, A_ref, n, info_ref)
                    if (uplos(ku) == 'U') then
                        do jg = 1, n
                            do ig = jg + 1, n
                                A_ref(ig, jg) = (0.0_ep, 0.0_ep)
                                A_got(ig, jg) = (0.0_ep, 0.0_ep)
                            end do
                        end do
                    else
                        do jg = 1, n
                            do ig = 1, jg - 1
                                A_ref(ig, jg) = (0.0_ep, 0.0_ep)
                                A_got(ig, jg) = (0.0_ep, 0.0_ep)
                            end do
                        end do
                    end if
                    err = max_rel_err_mat_z(A_got, A_ref)
                    tol = 64.0_ep * real(n, ep)**2 * target_eps
                    write(label, '(a,a1,a,a1,a,i0)') 'uplo=', uplos(ku), &
                        ',diag=', diags(kd), ',n=', n
                    call report_case(trim(label), err, tol)
                    deallocate(A_ref, A_got)
                end if
                deallocate(A_loc, A_glob)
            end do
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrtri
