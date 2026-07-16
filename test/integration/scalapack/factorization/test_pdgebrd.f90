program test_pdgebrd
    ! Bidiagonal reduction A -> Q*B*P^T. Compare the gathered post-call
    ! A (which encodes B + Householder reflectors) against dgebrd.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgebrd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgebrd
    implicit none

    integer, parameter :: ms(*) = [48, 64, 96]
    integer, parameter :: ns(*) = [32, 64, 64]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info, info_ref, lwork, mn
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: d(:), e(:), tauq(:), taup(:), work(:)
    real(ep), allocatable :: d_ref(:), e_ref(:), tauq_ref(:), taup_ref(:), work_ref(:)
    real(ep) :: err, tol

    call grid_init()
    call report_init('pdgebrd', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i); mn = min(m, n)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 13201 + 31*i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(d(max(1, locn_a)), e(max(1, locm_a)), &
                 tauq(max(1, locn_a)), taup(max(1, locm_a)), work(1))
        call target_pdgebrd(m, n, A_loc, 1, 1, desca, d, e, tauq, taup, &
                            work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdgebrd(m, n, A_loc, 1, 1, desca, d, e, tauq, taup, &
                            work, lwork, info)
        call gather_matrix(m, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n), d_ref(mn), e_ref(mn), &
                     tauq_ref(mn), taup_ref(mn), &
                     work_ref(max(1, (m + n) * 64)))
            A_ref = A_glob
            call dgebrd(m, n, A_ref, m, d_ref, e_ref, tauq_ref, taup_ref, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat(A_got, A_ref)
            tol = 64.0_ep * real(max(m, n), ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
                call report_case(trim(label), err, tol)
            end block
            deallocate(A_ref, d_ref, e_ref, tauq_ref, taup_ref, work_ref, A_got)
        end if
        deallocate(A_loc, A_glob, d, e, tauq, taup, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgebrd
