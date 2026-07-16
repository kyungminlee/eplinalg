program test_pzpoequ
    ! Z-mirror of pdpoequ. Hermitian-PD; scaling factors are real.
    use prec_kinds,       only: ep
    use compare,          only: rel_err_scalar, max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zpoequ
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_vector, &
                                gather_vector_row
    use target_scalapack, only: target_name, target_eps, target_pzpoequ
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), M_glob(:,:), dummy(:,:)
    real(ep), allocatable :: SR_loc(:), SC_loc(:), SR_got(:), SC_got(:), S_ref(:)
    real(ep) :: scond, amax, scond_ref, amax_ref, err, tol
    character(len=48) :: label
    integer :: ig, jg, owner_r, owner_c, il, jl, k

    call grid_init()
    call report_init('pzpoequ', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, dummy, M_glob, seed = 22201 + 31*i)
        deallocate(dummy)
        allocate(A_glob(n, n))
        A_glob = matmul(conjg(transpose(M_glob)), M_glob)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + cmplx(real(n, ep), 0.0_ep, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        allocate(A_loc(max(1, locm_a), max(1, locn_a)))
        A_loc = (0.0_ep, 0.0_ep)
        do jg = 1, n
            call g2l(jg, nb, my_npcol, owner_c, jl)
            if (owner_c /= my_col) cycle
            do ig = 1, n
                call g2l(ig, mb, my_nprow, owner_r, il)
                if (owner_r == my_row) A_loc(il, jl) = A_glob(ig, jg)
            end do
        end do
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(SR_loc(max(1, locm_a)), SC_loc(max(1, locn_a)))
        call target_pzpoequ(n, A_loc, 1, 1, desca, SR_loc, SC_loc, &
                            scond, amax, info)
        call gather_vector(n, mb, SR_loc, SR_got)
        call gather_vector_row(n, nb, SC_loc, SC_got)

        if (my_rank == 0) then
            allocate(S_ref(n))
            call zpoequ(n, A_glob, n, S_ref, scond_ref, amax_ref, info_ref)
            err = max(max_rel_err_vec(SR_got, S_ref), &
                      max_rel_err_vec(SC_got, S_ref), &
                      rel_err_scalar(scond, scond_ref), &
                      rel_err_scalar(amax, amax_ref))
            tol = 32.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(S_ref, SR_got, SC_got)
        end if
        deallocate(A_loc, A_glob, M_glob, SR_loc, SC_loc)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzpoequ
