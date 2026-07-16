program test_pzlanhs
    ! Z-mirror of test_pdlanhs.
    use prec_kinds,       only: ep
    use compare,          only: rel_err_scalar
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zlanhs
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzlanhs
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: norms(*) = [character(len=1) :: '1', 'I', 'F', 'M']
    integer :: i, j, n, info
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    real(ep),    allocatable :: work(:), work_ref(:)
    real(ep) :: got, refv, err, tol
    character(len=48) :: label
    integer :: ig, jg, owner_r, owner_c, il, jl

    call grid_init()
    call report_init('pzlanhs', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 9921 + 31*i)

        do jg = 1, n
            do ig = jg + 2, n
                A_glob(ig, jg) = (0.0_ep, 0.0_ep)
                call g2l(ig, mb, my_nprow, owner_r, il)
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_r == my_row .and. owner_c == my_col) &
                    A_loc(il, jl) = (0.0_ep, 0.0_ep)
            end do
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(work(max(1, max(locm_a, locn_a))))
        allocate(work_ref(max(1, n)))

        do j = 1, size(norms)
            got = target_pzlanhs(norms(j), n, A_loc, 1, 1, desca, work)
            if (my_rank == 0) then
                refv = zlanhs(norms(j), n, A_glob, n, work_ref)
                err = rel_err_scalar(got, refv)
                tol = 32.0_ep * real(n, ep) * target_eps
                write(label, '(a,a1,a,i0)') 'norm=', norms(j), ',n=', n
                call report_case(trim(label), err, tol)
            end if
        end do

        deallocate(A_loc, A_glob, work, work_ref)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzlanhs
