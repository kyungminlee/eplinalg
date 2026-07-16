program test_pzgetri
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgetrf, zgetri
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, &
                                target_pzgetrf, target_pzgetri
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, liwork, k
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    complex(ep), allocatable :: work(:), work_ref(:)
    integer,  allocatable :: ipiv_got(:), ipiv_ref(:), iwork(:)
    real(ep) :: err, tol
    complex(ep) :: boost

    call grid_init()
    call report_init('pzgetri', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 12501 + 31*i)

        boost = cmplx(real(n, ep), 0.0_ep, kind=ep)
        do k = 1, n
            A_glob(k, k) = A_glob(k, k) + boost
            call g2l(k, mb, my_nprow, owner_r, il)
            call g2l(k, nb, my_npcol, owner_c, jl)
            if (owner_r == my_row .and. owner_c == my_col) then
                A_loc(il, jl) = A_loc(il, jl) + boost
            end if
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(ipiv_got(locm_a + mb))
        call target_pzgetrf(n, n, A_loc, 1, 1, desca, ipiv_got, info)

        liwork = max(1, locn_a + mb)
        allocate(work(1), iwork(liwork))
        call target_pzgetri(n, A_loc, 1, 1, desca, ipiv_got, work, -1, &
                            iwork, -1, info)
        lwork  = max(1, int(real(work(1), ep)))
        liwork = max(liwork, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        call target_pzgetri(n, A_loc, 1, 1, desca, ipiv_got, work, lwork, &
                            iwork, liwork, info)
        call gather_matrix_z(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), ipiv_ref(n), work_ref(max(1, n * 64)))
            A_ref = A_glob
            call zgetrf(n, n, A_ref, n, ipiv_ref, info_ref)
            call zgetri(n, A_ref, n, ipiv_ref, work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(A_got, A_ref)
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
end program test_pzgetri
