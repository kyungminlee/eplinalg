program test_pdggrqf
    ! Generalized RQ factorization: A=R*Q (m×n), B=Z*T*Q (p×n).
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dggrqf
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local
    use pblas_distrib,     only: gen_distrib_matrix, gather_matrix
    use target_scalapack,  only: target_name, target_eps, target_pdggrqf
    implicit none

    integer, parameter :: m_v(*) = [32, 48, 48]
    integer, parameter :: p_v(*) = [40, 56, 64]
    integer, parameter :: n_v(*) = [48, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, p, n, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:), A_ref(:,:)
    real(ep), allocatable :: B_loc(:,:), B_glob(:,:), B_got(:,:), B_ref(:,:)
    real(ep), allocatable :: taua(:), taub(:), work(:)
    real(ep), allocatable :: taua_ref(:), taub_ref(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdggrqf', target_name, my_rank)

    do i = 1, size(m_v)
        m = m_v(i); p = p_v(i); n = n_v(i)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 27201 + 31*i)
        call gen_distrib_matrix(p, n, mb, nb, B_loc, B_glob, seed = 27211 + 31*i)

        locm_a = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(p, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(n, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
        call descinit_local(desca, m, n, mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, p, n, mb, nb, 0, 0, my_context, lld_b, info)

        ! pdggrqf TAUA is row-distributed of size LOCr(IA+M-1); TAUB is
        ! column-distributed of size LOCc(JB+MIN(P,N)-1). Use locn_b
        ! (= LOCc(N) >= LOCc(MIN(P,N))) for the column axis.
        allocate(taua(max(1, locm_a)), taub(max(1, locn_b)), work(1))
        call target_pdggrqf(m, p, n, A_loc, 1, 1, desca, taua, &
                            B_loc, 1, 1, descb, taub, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdggrqf(m, p, n, A_loc, 1, 1, desca, taua, &
                            B_loc, 1, 1, descb, taub, work, lwork, info)
        call gather_matrix(m, n, mb, nb, A_loc, A_got)
        call gather_matrix(p, n, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n), B_ref(p, n), taua_ref(min(m, n)), &
                     taub_ref(min(p, n)), &
                     work_ref(max(1, max(m, n, p) * 64)))
            A_ref = A_glob; B_ref = B_glob
            call dggrqf(m, p, n, A_ref, m, taua_ref, B_ref, p, taub_ref, &
                        work_ref, size(work_ref), info_ref)
            err = max(max_rel_err_mat(A_got, A_ref), &
                      max_rel_err_mat(B_got, B_ref))
            tol = 64.0_ep * real(max(m, n, p), ep)**2 * target_eps
            write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',p=', p, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, taua_ref, taub_ref, work_ref, A_got, B_got)
        end if
        deallocate(A_loc, A_glob, B_loc, B_glob, taua, taub, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdggrqf
