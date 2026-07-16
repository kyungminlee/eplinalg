program test_pztrsm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: ztrsm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local, g2l
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pztrsm
    implicit none

    integer, parameter :: ms(*) = [32, 64, 128]
    integer, parameter :: ns(*) = [40, 60, 100]
    ! Cover the 5 representative SIDE/UPLO/TRANSA/DIAG combos plus a
    ! TRANSA='C' case and an alpha != 1 scaling case.
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'L', 'R', 'R', 'L', 'L']
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'U', 'L', 'L', 'U']
    character(len=1), parameter :: transas(*) = ['N', 'N', 'T', 'N', 'N', 'C', 'N']
    character(len=1), parameter :: diags(*)   = ['N', 'N', 'N', 'N', 'U', 'N', 'N']
    real(ep), parameter :: alpha_re(*) = [1.0_ep, 1.0_ep, 1.0_ep, 1.0_ep, 1.0_ep, 1.0_ep, 0.6_ep]
    real(ep), parameter :: alpha_im(*) = [0.0_ep, 0.0_ep, 0.0_ep, 0.0_ep, 0.0_ep, 0.0_ep, 0.2_ep]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, j, ic, m, n, info
    integer :: ka
    integer :: locm_a, locn_a, locm_b, locn_b, lld_a, lld_b
    integer :: desca(9), descb(9)
    integer :: owner_r, owner_c, il, jl
    character(len=1) :: side, uplo, transa, diag
    complex(ep) :: bump, alpha
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pztrsm', target_name, my_rank)

    do ic = 1, size(sides)
        side = sides(ic); uplo = uplos(ic)
        transa = transas(ic); diag = diags(ic)
        alpha = cmplx(alpha_re(ic), alpha_im(ic), ep)
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (side == 'L') then
                ka = m
            else
                ka = n
            end if

            call gen_distrib_matrix_z(ka, ka, mb, mb, A_loc, A_glob, &
                                      seed = 14701 + 59 * i + 311 * ic)
            call gen_distrib_matrix_z(m, n, mb, nb, B_loc, B0, &
                                      seed = 14711 + 59 * i + 311 * ic)

            bump = cmplx(real(ka, ep), 0.0_ep, ep)
            do j = 1, ka
                A_glob(j, j) = A_glob(j, j) + bump
                call g2l(j, mb, my_nprow, owner_r, il)
                call g2l(j, mb, my_npcol, owner_c, jl)
                if (owner_r == my_row .and. owner_c == my_col) then
                    A_loc(il, jl) = A_loc(il, jl) + bump
                end if
            end do

            locm_a = numroc_local(ka, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(ka, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_b = numroc_local(m, mb, my_row, 0, my_nprow)
            locn_b = numroc_local(n, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)

            call descinit_local(desca, ka, ka, mb, mb, 0, 0, my_context, lld_a, info)
            call descinit_local(descb, m, n, mb, nb, 0, 0, my_context, lld_b, info)

            call target_pztrsm(side, uplo, transa, diag, m, n, alpha, &
                               A_loc, 1, 1, desca, B_loc, 1, 1, descb)
            call gather_matrix_z(m, n, mb, nb, B_loc, B_got)

            if (my_rank == 0) then
                allocate(B_ref(m, n))
                B_ref = B0
                call ztrsm(side, uplo, transa, diag, m, n, alpha, &
                           A_glob, ka, B_ref, m)
                err = max_rel_err_mat_z(B_got, B_ref)
                tol = 64.0_ep * 8.0_ep * real(ka, ep) * target_eps
                write(label, '(a,a,a,a,a,a,a,a,a,i0,a,i0)') &
                    'side=', side, ',uplo=', uplo, ',ta=', transa, ',d=', diag, &
                    ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(B_ref, B_got)
            end if
            deallocate(A_loc, B_loc, A_glob, B0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrsm
