program test_pdlantr
    use prec_kinds,       only: ep
    use compare,          only: rel_err_scalar
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dlantr
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix
    use target_scalapack, only: target_name, target_eps, target_pdlantr
    implicit none

    integer, parameter :: ms(*) = [40, 64, 96]
    integer, parameter :: ns(*) = [32, 64, 80]
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: norms(*) = [character(len=1) :: '1', 'I', 'F', 'M']
    character(len=1), parameter :: uplos(*) = [character(len=1) :: 'U', 'L']
    character(len=1), parameter :: diags(*) = [character(len=1) :: 'N', 'U']
    integer :: i, j, k, l, m, n, info
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    real(ep) :: got, refv, err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdlantr', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 15601 + 31*i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(work(max(1, max(m, n))), work_ref(max(1, max(m, n))))

        do l = 1, size(diags)
            do k = 1, size(uplos)
                do j = 1, size(norms)
                    got = target_pdlantr(norms(j), uplos(k), diags(l), m, n, &
                                         A_loc, 1, 1, desca, work)
                    if (my_rank == 0) then
                        refv = dlantr(norms(j), uplos(k), diags(l), m, n, &
                                      A_glob, m, work_ref)
                        err = rel_err_scalar(got, refv)
                        tol = 32.0_ep * real(max(m, n), ep) * target_eps
                        write(label, '(a,a1,a,a1,a,a1,a,i0,a,i0)') &
                            'norm=', norms(j), ',uplo=', uplos(k), &
                            ',diag=', diags(l), ',m=', m, ',n=', n
                        call report_case(trim(label), err, tol)
                    end if
                end do
            end do
        end do

        deallocate(A_loc, A_glob, work, work_ref)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdlantr
