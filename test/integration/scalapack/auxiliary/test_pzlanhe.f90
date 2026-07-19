program test_pzlanhe
    use prec_kinds,       only: ep
    use compare,          only: rel_err_scalar
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zlanhe
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix_z, scatter_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzlanhe
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: norms(*) = [character(len=1) :: '1', 'I', 'F', 'M']
    character(len=1), parameter :: uplos(*) = [character(len=1) :: 'U', 'L']
    integer :: i, j, k, n, info
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_herm(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    real(ep) :: got, refv, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzlanhe', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 15801 + 31*i)
        allocate(A_herm(n, n))
        A_herm = 0.5_ep * (A_glob + conjg(transpose(A_glob)))

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call scatter_matrix_z(n, n, mb, nb, A_herm, A_loc)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(work(max(1024, 8 * n)), work_ref(max(1024, 8 * n)))

        do k = 1, size(uplos)
            do j = 1, size(norms)
                got = target_pzlanhe(norms(j), uplos(k), n, A_loc, 1, 1, desca, work)
                if (my_rank == 0) then
                    refv = zlanhe(norms(j), uplos(k), n, A_herm, n, work_ref)
                    err = rel_err_scalar(got, refv)
                    tol = 32.0_ep * real(n, ep) * target_eps
                    write(label, '(a,a1,a,a1,a,i0)') 'norm=', norms(j), &
                        ',uplo=', uplos(k), ',n=', n
                    call report_case(trim(label), err, tol)
                end if
            end do
        end do

        deallocate(A_loc, A_glob, A_herm, work, work_ref)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzlanhe
