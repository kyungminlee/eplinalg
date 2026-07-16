program test_pztrsv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: ztrsv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local, g2l
    use pblas_distrib, only: gen_distrib_matrix_z, gen_distrib_vector_z, &
                             gather_vector_z
    use target_pblas,  only: target_name, target_eps, target_pztrsv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    character(len=1), parameter :: uplos(*)   = ['U', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'T', 'C']
    character(len=1), parameter :: diags(*)   = ['N', 'U']
    integer, parameter :: mb = 8
    integer :: i, iu, it, id, j, n, info
    integer :: combo
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9)
    integer :: owner_r, owner_c, il, jl
    character(len=1) :: uplo, trans, diag
    complex(ep) :: bump
    complex(ep), allocatable :: A_loc(:,:), x_loc(:)
    complex(ep), allocatable :: A_glob(:,:), x_glob(:), x_got(:), x_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pztrsv', target_name, my_rank)

    do iu = 1, size(uplos)
        do it = 1, size(transes)
            do id = 1, size(diags)
                uplo = uplos(iu); trans = transes(it); diag = diags(id)
                combo = ((iu - 1) * size(transes) + (it - 1)) * size(diags) + id
                do i = 1, size(ns)
                    n = ns(i)
                    call gen_distrib_matrix_z(n, n, mb, mb, A_loc, A_glob, &
                                              seed = 6001 + 19 * i + 211 * combo)
                    call gen_distrib_vector_z(n, mb, x_loc, x_glob, &
                                              seed = 6011 + 19 * i + 211 * combo)

                    bump = cmplx(real(n, ep), 0.0_ep, ep)
                    do j = 1, n
                        A_glob(j, j) = A_glob(j, j) + bump
                        call g2l(j, mb, my_nprow, owner_r, il)
                        call g2l(j, mb, my_npcol, owner_c, jl)
                        if (owner_r == my_row .and. owner_c == my_col) then
                            A_loc(il, jl) = A_loc(il, jl) + bump
                        end if
                    end do

                    locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
                    locn_a = numroc_local(n, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
                    locn_x = numroc_local(n, mb, my_row, 0, my_nprow); lld_x = max(1, locn_x)

                    call descinit_local(desca, n, n, mb, mb, 0, 0, my_context, lld_a, info)
                    call descinit_local(descx, n, 1, mb, 1, 0, 0, my_context, lld_x, info)

                    call target_pztrsv(uplo, trans, diag, n, A_loc, 1, 1, desca, &
                                       x_loc, 1, 1, descx, 1)
                    call gather_vector_z(n, mb, x_loc, x_got)

                    if (my_rank == 0) then
                        allocate(x_ref(n))
                        x_ref = x_glob
                        call ztrsv(uplo, trans, diag, n, A_glob, n, x_ref, 1)
                        err = max_rel_err_vec_z(x_got, x_ref)
                        tol = 64.0_ep * 8.0_ep * real(n, ep) * target_eps
                        write(label, '(a,a,a,a,a,a,a,i0)') &
                            'uplo=', uplo, ',trans=', trans, ',diag=', diag, ',n=', n
                        call report_case(trim(label), err, tol)
                        deallocate(x_ref, x_got)
                    end if
                    deallocate(A_loc, x_loc, A_glob, x_glob)
                end do
            end do
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrsv
