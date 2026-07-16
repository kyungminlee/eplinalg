program test_pzhemm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zhemm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzhemm
    implicit none

    integer, parameter :: ms(*) = [32, 64, 128]
    integer, parameter :: ns(*) = [32, 48, 96]
    character(len=1), parameter :: sides(*) = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U', 'L']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, is, m, n, info
    integer :: ka
    integer :: locm_a, locn_a, locm_b, locn_b, locm_c, locn_c
    integer :: lld_a, lld_b, lld_c
    integer :: desca(9), descb(9), descc(9)
    character(len=1) :: side, uplo
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), &
                                C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzhemm', target_name, my_rank)

    alpha = cmplx(0.7_ep, 0.2_ep, ep); beta = cmplx(0.3_ep, -0.1_ep, ep)
    do is = 1, size(sides)
        side = sides(is); uplo = uplos(is)
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (side == 'L') then
                ka = m
            else
                ka = n
            end if

            call gen_distrib_matrix_z(ka, ka, mb, mb, A_loc, A_glob, &
                                      seed = 14101 + 47 * i + 211 * is)
            call gen_distrib_matrix_z(m, n, mb, nb, B_loc, B_glob, &
                                      seed = 14111 + 47 * i + 211 * is)
            call gen_distrib_matrix_z(m, n, mb, nb, C_loc, C0, &
                                      seed = 14121 + 47 * i + 211 * is)

            locm_a = numroc_local(ka, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(ka, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_b = numroc_local(m, mb, my_row, 0, my_nprow)
            locn_b = numroc_local(n, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
            locm_c = locm_b; locn_c = locn_b; lld_c = lld_b

            call descinit_local(desca, ka, ka, mb, mb, 0, 0, my_context, lld_a, info)
            call descinit_local(descb, m, n, mb, nb, 0, 0, my_context, lld_b, info)
            call descinit_local(descc, m, n, mb, nb, 0, 0, my_context, lld_c, info)

            call target_pzhemm(side, uplo, m, n, alpha, &
                               A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                               beta,  C_loc, 1, 1, descc)
            call gather_matrix_z(m, n, mb, nb, C_loc, C_got)

            if (my_rank == 0) then
                allocate(C_ref(m, n))
                C_ref = C0
                call zhemm(side, uplo, m, n, alpha, A_glob, ka, B_glob, m, &
                           beta, C_ref, m)
                err = max_rel_err_mat_z(C_got, C_ref)
                tol = 64.0_ep * 8.0_ep * real(ka, ep) * target_eps
                write(label, '(a,a,a,a,a,i0,a,i0)') &
                    'side=', side, ',uplo=', uplo, ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(C_ref, C_got)
            end if
            deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzhemm
