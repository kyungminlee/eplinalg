program test_pdorg2r
    ! Build Q from a QR factorization computed by pdgeqrf, then compare
    ! against the reference Q from dgeqrf+dorgqr. Verifies Q itself, not
    ! just the orthogonality residual.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgeqrf, dorgqr
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdgeqrf, target_pdorg2r
    implicit none

    integer, parameter :: ms(*) = [48, 80, 96]
    integer, parameter :: ns(*) = [32, 48, 64]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, k, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), Q_got(:,:), A_ref(:,:)
    real(ep), allocatable :: tau_got(:), tau_ref(:), work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdorg2r', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i); k = n
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 9501 + 31*i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)

        ! QR factorize first.
        allocate(tau_got(max(1, locn_a)), work(1))
        call target_pdgeqrf(m, n, A_loc, 1, 1, desca, tau_got, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdgeqrf(m, n, A_loc, 1, 1, desca, tau_got, work, lwork, info)
        deallocate(work)

        ! Generate Q.
        allocate(work(1))
        call target_pdorg2r(m, n, k, A_loc, 1, 1, desca, tau_got, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdorg2r(m, n, k, A_loc, 1, 1, desca, tau_got, work, lwork, info)
        call gather_matrix(m, n, mb, nb, A_loc, Q_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n), tau_ref(min(m, n)), work_ref(max(1, n * 64)))
            A_ref = A_glob
            call dgeqrf(m, n, A_ref, m, tau_ref, work_ref, size(work_ref), info_ref)
            call dorgqr(m, n, k, A_ref, m, tau_ref, work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat(Q_got, A_ref)
            tol = 64.0_ep * real(max(m, n), ep)**2 * target_eps
            write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(A_ref, tau_ref, work_ref, Q_got)
        end if
        deallocate(A_loc, A_glob, tau_got, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdorg2r
