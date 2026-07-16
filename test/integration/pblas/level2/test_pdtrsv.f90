program test_pdtrsv
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_vec
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dtrsv
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local, g2l
    use pblas_distrib, only: gen_distrib_matrix, gen_distrib_vector, &
                             gather_vector
    use target_pblas,  only: target_name, target_eps, target_pdtrsv
    implicit none

    integer, parameter :: ns(*) = [32, 80, 160]
    character(len=1), parameter :: uplos(*)   = ['U', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'T']
    character(len=1), parameter :: diags(*)   = ['N', 'U']
    integer, parameter :: mb = 8
    integer :: i, iu, it, id, j, n, info
    integer :: combo
    integer :: locm_a, locn_a, locn_x, lld_a, lld_x
    integer :: desca(9), descx(9)
    integer :: owner_r, owner_c, il, jl
    character(len=1) :: uplo, trans, diag
    real(ep) :: bump
    real(ep), allocatable :: A_loc(:,:), x_loc(:)
    real(ep), allocatable :: A_glob(:,:), x_glob(:), x_got(:), x_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdtrsv', target_name, my_rank)

    do iu = 1, size(uplos)
        do it = 1, size(transes)
            do id = 1, size(diags)
                uplo = uplos(iu); trans = transes(it); diag = diags(id)
                combo = ((iu - 1) * size(transes) + (it - 1)) * size(diags) + id
                do i = 1, size(ns)
                    n = ns(i)
                    call gen_distrib_matrix(n, n, mb, mb, A_loc, A_glob, &
                                            seed = 3101 + 19 * i + 211 * combo)
                    call gen_distrib_vector(n, mb, x_loc, x_glob, &
                                            seed = 3111 + 19 * i + 211 * combo)

                    ! Diagonal bump for diagonal-dominance — relevant only
                    ! for DIAG='N' (DIAG='U' ignores diagonal entirely),
                    ! but harmless to apply unconditionally. Works for
                    ! both UPLO='U' and UPLO='L' since the read triangle
                    ! always includes A(j,j).
                    bump = real(n, ep)
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

                    call target_pdtrsv(uplo, trans, diag, n, A_loc, 1, 1, desca, &
                                       x_loc, 1, 1, descx, 1)
                    call gather_vector(n, mb, x_loc, x_got)

                    if (my_rank == 0) then
                        allocate(x_ref(n))
                        x_ref = x_glob
                        call dtrsv(uplo, trans, diag, n, A_glob, n, x_ref, 1)
                        err = max_rel_err_vec(x_got, x_ref)
                        tol = 64.0_ep * real(n, ep) * target_eps
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
end program test_pdtrsv
