program test_pdtradd
    ! pdtradd: trapezoidal C := beta*C + alpha*op(A). Reference is
    ! element-wise inside the named triangle of C.
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix, gather_matrix
    use target_pblas,  only: target_name, target_eps, target_pdtradd
    implicit none

    integer, parameter :: cases = 4
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'L']
    character(len=1), parameter :: transes(*) = ['N', 'N', 'T', 'T']
    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, ic, n, info, ig, jg
    integer :: locm_a, locn_a, locm_c, locn_c, lld_a, lld_c
    integer :: desca(9), descc(9)
    real(ep), allocatable :: A_loc(:,:), C_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    real(ep) :: alpha, beta, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdtradd', target_name, my_rank)

    alpha = 0.7_ep; beta = 0.3_ep
    do ic = 1, cases
        do i = 1, size(ns)
            n = ns(i)
            call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, &
                                    seed = 17201 + 29*i + 211*ic)
            call gen_distrib_matrix(n, n, mb, nb, C_loc, C0, &
                                    seed = 17211 + 29*i + 211*ic)

            locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_c = locm_a; locn_c = locn_a; lld_c = lld_a
            call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)
            call descinit_local(descc, n, n, mb, nb, 0, 0, my_context, lld_c, info)

            call target_pdtradd(uplos(ic), transes(ic), n, n, alpha, &
                                A_loc, 1, 1, desca, beta, C_loc, 1, 1, descc)
            call gather_matrix(n, n, mb, nb, C_loc, C_got)

            if (my_rank == 0) then
                allocate(C_ref(n, n))
                C_ref = C0
                if (transes(ic) == 'N') then
                    C_ref = beta * C_ref + alpha * A_glob
                else
                    C_ref = beta * C_ref + alpha * transpose(A_glob)
                end if
                ! Mask off the off-triangle so comparison only sees the
                ! piece both target and reference write to.
                if (uplos(ic) == 'U') then
                    do jg = 1, n
                        do ig = jg + 1, n
                            C_ref(ig, jg) = 0.0_ep
                            C_got(ig, jg) = 0.0_ep
                        end do
                    end do
                else
                    do jg = 1, n
                        do ig = 1, jg - 1
                            C_ref(ig, jg) = 0.0_ep
                            C_got(ig, jg) = 0.0_ep
                        end do
                    end do
                end if
                err = max_rel_err_mat(C_got, C_ref)
                tol = 16.0_ep * real(n, ep) * target_eps
                write(label, '(a,a,a,a,a,i0)') 'uplo=', uplos(ic), &
                    ',trans=', transes(ic), ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(C_ref, C_got)
            end if
            deallocate(A_loc, C_loc, A_glob, C0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pdtradd
