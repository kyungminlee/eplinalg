program test_pdgetri
    ! LU factorize then invert; compare against dgetrf+dgetri.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgetrf, dgetri
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdgetrf, target_pdgetri
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: ipiv_got(:), ipiv_ref(:), iwork(:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pdgetri', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 12201 + 31*i)

        ! Diagonal boost for conditioning.
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + real(n, ep)
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + real(n, ep)
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(ipiv_got(locm_a + mb))
        call target_pdgetrf(n, n, A_loc, 1, 1, desca, ipiv_got, info)

        ! Workspace query.
        liwork = max(1, locn_a + mb)
        allocate(work(1), iwork(liwork))
        call target_pdgetri(n, A_loc, 1, 1, desca, ipiv_got, work, -1, &
                            iwork, -1, info)
        lwork = max(1, int(work(1)))
        liwork = max(liwork, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pdgetri(n, A_loc, 1, 1, desca, ipiv_got, work, lwork, &
                            iwork, liwork, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), ipiv_ref(n), work_ref(max(1, n * 64)))
            A_ref = A_glob
            call dgetrf(n, n, A_ref, n, ipiv_ref, info_ref)
            call dgetri(n, A_ref, n, ipiv_ref, work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,i0)') 'n=', n
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, ipiv_ref, work_ref, A_got)
        end if
        deallocate(A_loc, A_glob, ipiv_got, work, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgetri
