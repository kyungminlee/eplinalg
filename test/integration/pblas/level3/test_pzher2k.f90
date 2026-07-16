! NOTE: this test currently *fails* on every staging build because the
! migrated PB_Cconjg / PB_Ctzher2k stack-allocate ``cmplx16 Calph16``
! (16 bytes, 2 × ``double``) and the C migrator does not rewrite
! ``cmplx16`` to ``cmplxQ`` in local-variable declarations. The
! complex alpha conjugate written by PB_Cconjg into Calph16 is
! truncated, and the subsequent gemm call reads garbage for
! conjg(alpha). See ``tests/pblas/TODO.md`` for the full diagnosis
! and the recommended fix in ``codegen/migrator/c_migrator.py``. The test
! is left in place so the failure is visible until the migrator fix
! lands; a real-only alpha would mask the bug and is not what we want.
!
! Additional UPLO/TRANS combinations have been added alongside the
! original (UPLO='U', TRANS='N') case for path coverage. They route
! through the same PB_Ctzher2k/PB_Cconjg stack-allocation site, so they
! also fail on the kind16 target until the migrator fix lands. They
! are intentionally left in place so a fix can be validated across all
! four PB_Cptzher2k{LU,LL,UN,UC,LN,LC} paths at once.

program test_pzher2k
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report,   only: report_init, report_case, report_finalize
    use pblas_ref_quad_blas, only: zher2k
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzher2k
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: ks(*) = [20, 40, 72]
    ! TRANS='T' is invalid for her2k; 'C' is the conjugate-transpose path.
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'C', 'C']
    integer, parameter :: mb = 8
    integer :: i, ic, n, k, info
    integer :: ar, ac
    integer :: locm_a, locn_a, locm_b, locn_b, locm_c, locn_c
    integer :: lld_a, lld_b, lld_c
    integer :: desca(9), descb(9), descc(9)
    character(len=1) :: uplo, trans
    complex(ep), allocatable :: A_loc(:,:), B_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), B_glob(:,:), C0(:,:), &
                                C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha
    real(ep) :: beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzher2k', target_name, my_rank)

    alpha = cmplx(0.5_ep, 0.2_ep, ep); beta = 0.25_ep
    do ic = 1, size(uplos)
        uplo = uplos(ic); trans = transes(ic)
        do i = 1, size(ns)
            n = ns(i); k = ks(i)
            if (trans == 'N') then
                ar = n; ac = k
            else
                ar = k; ac = n
            end if

            call gen_distrib_matrix_z(ar, ac, mb, mb, A_loc, A_glob, &
                                      seed = 14501 + 53 * i + 211 * ic)
            call gen_distrib_matrix_z(ar, ac, mb, mb, B_loc, B_glob, &
                                      seed = 14511 + 53 * i + 211 * ic)
            call gen_distrib_matrix_z(n, n, mb, mb, C_loc, C0, &
                                      seed = 14521 + 53 * i + 211 * ic)

            locm_a = numroc_local(ar, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(ac, mb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_b = locm_a; locn_b = locn_a; lld_b = lld_a
            locm_c = numroc_local(n, mb, my_row, 0, my_nprow)
            locn_c = numroc_local(n, mb, my_col, 0, my_npcol); lld_c = max(1, locm_c)

            call descinit_local(desca, ar, ac, mb, mb, 0, 0, my_context, lld_a, info)
            call descinit_local(descb, ar, ac, mb, mb, 0, 0, my_context, lld_b, info)
            call descinit_local(descc, n, n, mb, mb, 0, 0, my_context, lld_c, info)

            call target_pzher2k(uplo, trans, n, k, alpha, &
                                A_loc, 1, 1, desca, B_loc, 1, 1, descb, &
                                beta, C_loc, 1, 1, descc)
            call gather_matrix_z(n, n, mb, mb, C_loc, C_got)

            if (my_rank == 0) then
                allocate(C_ref(n, n))
                C_ref = C0
                call zher2k(uplo, trans, n, k, alpha, A_glob, ar, B_glob, ar, &
                            beta, C_ref, n)
                err = max_rel_err_mat_z(C_got, C_ref)
                tol = 64.0_ep * 8.0_ep * real(k, ep) * target_eps
                write(label, '(a,a,a,a,a,i0,a,i0)') &
                    'uplo=', uplo, ',trans=', trans, ',n=', n, ',k=', k
                call report_case(trim(label), err, tol)
                deallocate(C_ref, C_got)
            end if
            deallocate(A_loc, B_loc, C_loc, A_glob, B_glob, C0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzher2k
