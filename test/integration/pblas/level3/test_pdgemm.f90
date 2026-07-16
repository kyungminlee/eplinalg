program test_pdgemm
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: dgemm
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdgemm
    implicit none

    ! ms/ns/ks include one shape (m=70, n=55, k=43) that is *not* a
    ! multiple of mb=nb=8 — the partial-block edge is the most likely
    ! place for a numroc/g2l off-by-one to surface, and every
    ! aligned-only shape would miss it.
    integer, parameter :: ms(*) = [32, 80, 160, 70, 32, 32, 64, 64]
    integer, parameter :: ns(*) = [40, 60, 120, 55, 40, 40, 48, 48]
    integer, parameter :: ks(*) = [24, 48, 100, 43, 24, 24, 40, 40]
    ! Per-shape (alpha, beta, transa, transb). Cases:
    !   1-4: TRANSA=N, TRANSB=N — main path with non-aligned shape (#4)
    !   5:   alpha=0 quick-return
    !   6:   beta=0 quick-return
    !   7:   TRANSA=T, TRANSB=N — exercises PB_CpgemmAB transpose paths
    !   8:   TRANSA=N, TRANSB=T
    real(ep), parameter :: alphas(*) = [0.7_ep, 0.7_ep, 0.7_ep, 0.7_ep, &
                                         0.0_ep, 0.5_ep, 0.6_ep, 0.6_ep]
    real(ep), parameter :: betas(*)  = [0.3_ep, 0.3_ep, 0.3_ep, 0.3_ep, &
                                         0.4_ep, 0.0_ep, 0.2_ep, 0.2_ep]
    character(len=1), parameter :: transas(*) = ['N','N','N','N','N','N','T','N']
    character(len=1), parameter :: transbs(*) = ['N','N','N','N','N','N','N','T']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, k, info
    integer :: ar, ac, br, bc
    integer :: locm_a, locn_a, locm_b, locn_b, locm_c, locn_c
    integer :: lld_a, lld_b, lld_c
    integer :: desca(9), descb(9), descc(9)
    character(len=1) :: transa, transb
    real(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pdgemm', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = ks(i)
        alpha = alphas(i); beta = betas(i)
        transa = transas(i); transb = transbs(i)
        ! Logical A is op(A) of size m×k; storage A is ar×ac.
        if (transa == 'N') then
            ar = m; ac = k
        else
            ar = k; ac = m
        end if
        if (transb == 'N') then
            br = k; bc = n
        else
            br = n; bc = k
        end if

        call gen_distrib_matrix(ar, ac, mb, nb, A_loc, A_glob, seed = 7801 + 23 * i)
        call gen_distrib_matrix(br, bc, mb, nb, B_loc, B_glob, seed = 7811 + 23 * i)
        call gen_distrib_matrix(m,  n,  mb, nb, C_loc, C0,     seed = 7821 + 23 * i)

        locm_a = numroc_local(ar, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(ac, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(br, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(bc, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        locm_c = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(n, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)

        call descinit_local(desca, ar, ac, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, br, bc, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descc, m,  n,  mb, nb, 0, 0, my_context, lld_c, info)

        call target_pdgemm(transa, transb, m, n, k, alpha, &
                           A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                           beta,  C_loc, 1, 1, descc)
        call gather_matrix(m, n, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(m, n))
            C_ref = C0
            call dgemm(transa, transb, m, n, k, alpha, A_glob, ar, B_glob, br, &
                       beta, C_ref, m)
            err = max_rel_err_mat(C_got, C_ref)
            tol = 32.0_ep * 2.0_ep * real(k, ep) * target_eps
            write(label, '(a,a,a,a,a,f3.1,a,f3.1,a,i0,a,i0,a,i0)') &
                'ta=', transa, ',tb=', transb, &
                ',a=', alpha, ',b=', beta, ',m=', m, ',n=', n, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgemm
