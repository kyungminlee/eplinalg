program test_pzgeadd
    ! pzgeadd: C := beta*C + alpha*op(A) for op in {N, T, C}.
    use prec_kinds,    only: ep
    use compare,       only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,    only: grid_init, grid_exit, my_rank, my_context, &
                             my_nprow, my_npcol, my_row, my_col, &
                             numroc_local, descinit_local
    use pblas_distrib, only: gen_distrib_matrix_z, gather_matrix_z
    use target_pblas,  only: target_name, target_eps, target_pzgeadd
    implicit none

    integer, parameter :: cases = 6
    character(len=1), parameter :: transes(*) = ['N', 'T', 'C', 'N', 'T', 'C']
    integer, parameter :: ms(*) = [32, 64, 96]
    integer, parameter :: ns(*) = [40, 48, 80]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, ic, m, n, ma, na_dim, info
    integer :: locm_a, locn_a, locm_c, locn_c, lld_a, lld_c
    integer :: desca(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), C_loc(:,:)
    complex(ep), allocatable :: A_glob(:,:), C0(:,:), C_ref(:,:), C_got(:,:)
    complex(ep) :: alpha, beta
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzgeadd', target_name, my_rank)

    alpha = cmplx(0.7_ep,  0.2_ep, kind=ep)
    beta  = cmplx(0.3_ep, -0.1_ep, kind=ep)
    do ic = 1, cases
        do i = 1, size(ms)
            m = ms(i); n = ns(i)
            if (transes(ic) == 'N') then
                ma = m; na_dim = n
            else
                ma = n; na_dim = m
            end if
            call gen_distrib_matrix_z(ma, na_dim, mb, nb, A_loc, A_glob, &
                                      seed = 17301 + 29*i + 211*ic)
            call gen_distrib_matrix_z(m, n, mb, nb, C_loc, C0, &
                                      seed = 17311 + 29*i + 211*ic)

            locm_a = numroc_local(ma, mb, my_row, 0, my_nprow)
            locn_a = numroc_local(na_dim, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
            locm_c = numroc_local(m, mb, my_row, 0, my_nprow)
            locn_c = numroc_local(n, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
            call descinit_local(desca, ma, na_dim, mb, nb, 0, 0, my_context, lld_a, info)
            call descinit_local(descc, m, n, mb, nb, 0, 0, my_context, lld_c, info)

            call target_pzgeadd(transes(ic), m, n, alpha, A_loc, 1, 1, desca, &
                                beta, C_loc, 1, 1, descc)
            call gather_matrix_z(m, n, mb, nb, C_loc, C_got)

            if (my_rank == 0) then
                allocate(C_ref(m, n))
                C_ref = C0
                select case (transes(ic))
                case ('N')
                    C_ref = beta * C_ref + alpha * A_glob
                case ('T')
                    C_ref = beta * C_ref + alpha * transpose(A_glob)
                case ('C')
                    C_ref = beta * C_ref + alpha * conjg(transpose(A_glob))
                end select
                err = max_rel_err_mat_z(C_got, C_ref)
                tol = 16.0_ep * real(max(m, n), ep) * target_eps
                write(label, '(a,a,a,i0,a,i0)') 'trans=', transes(ic), &
                    ',m=', m, ',n=', n
                call report_case(trim(label), err, tol)
                deallocate(C_ref, C_got)
            end if
            deallocate(A_loc, C_loc, A_glob, C0)
        end do
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgeadd
