program test_pdtrmm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dtrmm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdtrmm
    implicit none

    integer, parameter :: ms(*) = [32, 80, 128]
    integer, parameter :: ns(*) = [40, 60, 100]
    ! Cover the 5 representative SIDE/UPLO/TRANSA/DIAG combos called out
    ! in the test plan: L/U/N/N (original), L/L/N/N, L/U/T/N, R/U/N/N,
    ! R/L/N/U. Each case routes through a distinct PB_Cptrmm{...} path.
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'L', 'R', 'R']
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'U', 'L']
    character(len=1), parameter :: transas(*) = ['N', 'N', 'T', 'N', 'N']
    character(len=1), parameter :: diags(*)   = ['N', 'N', 'N', 'N', 'U']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, ic, m, n, info
    integer :: ka
    integer :: locm_a, locn_a, locm_b, locn_b, lld_a, lld_b
    integer :: desca(9), descb(9)
    character(len=1) :: side, uplo, transa, diag
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B0(:,:), B_ref(:,:), B_got(:,:)
    real(ep) :: alpha, err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pdtrmm', target_name, my_rank)

    alpha = 0.8_ep
    do ic = 1, size(sides)
        side = sides(ic); uplo = uplos(ic)
        transa = transas(ic); diag = diags(ic)
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (side == 'L') then
                ka = m
            else
                ka = n
            end if

            call gen_distrib_matrix(ka, ka, mb, mb, A_loc, A_glob, &
                                    seed = 10001 + 37 * i + 311 * ic)
            call gen_distrib_matrix(m, n, mb, nb, B_loc, B0, &
                                    seed = 10011 + 37 * i + 311 * ic)

            locm_a = numroc_local(ka, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(ka, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_b = numroc_local(m, mb, my_row, 0, my_nprow)
            locn_b = numroc_local(n, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)

            call descinit_local(desca, ka, ka, mb, mb, 0, 0, my_context, lld_a, info)
            call descinit_local(descb, m, n, mb, nb, 0, 0, my_context, lld_b, info)

            call target_pdtrmm(side, uplo, transa, diag, m, n, alpha, &
                               A_loc, 1, 1, desca, B_loc, 1, 1, descb)
            call gather_matrix(m, n, mb, nb, B_loc, B_got)

            if (my_rank == 0) then
                allocate(B_ref(m, n))
                B_ref = B0
                call dtrmm(side, uplo, transa, diag, m, n, alpha, &
                           A_glob, ka, B_ref, m)
                err = max_rel_err_mat(B_got, B_ref)
                tol = 32.0_ep * 2.0_ep * real(ka, ep) * target_eps
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
end program test_pdtrmm
