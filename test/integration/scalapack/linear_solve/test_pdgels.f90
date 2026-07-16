program test_pdgels
    ! Distributed least-squares: overdetermined m greater-than n with
    ! TRANS=N. Compare migrated result against serial dgels.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgels
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgels
    implicit none

    integer, parameter :: cases = 2
    character(len=1), parameter :: transes(*) = ['N', 'N']
    integer, parameter :: ms(*) = [64, 96]
    integer, parameter :: ns(*) = [32, 48]
    integer, parameter :: nrhs = 3
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, m, n, m_big, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a, locm_b, lld_b
    integer :: desca(9), descb(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), B_loc(:,:), B_glob(:,:)
    real(ep), allocatable :: B_got(:,:), A_ref(:,:), B_ref(:,:)
    real(ep), allocatable :: work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdgels', target_name, my_rank)

    do ic = 1, cases
        m = ms(ic); n = ns(ic); m_big = max(m, n)
        call gen_distrib_matrix(m,  n,    mb, nb, A_loc, A_glob, seed = 11001 + 31*ic)
        call gen_distrib_matrix(m_big, nrhs, mb, nb, B_loc, B_glob, seed = 11011 + 31*ic)

        locm_a = numroc_local(m,  mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n,  nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(m_big, mb, my_row, 0, my_nprow); lld_b = max(1, locm_b)
        call descinit_local(desca, m,  n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb, m_big, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(work(1))
        call target_pdgels(transes(ic), m, n, nrhs, A_loc, 1, 1, desca, &
                           B_loc, 1, 1, descb, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdgels(transes(ic), m, n, nrhs, A_loc, 1, 1, desca, &
                           B_loc, 1, 1, descb, work, lwork, info)
        call gather_matrix(m_big, nrhs, mb, nb, B_loc, B_got)

        if (my_rank == 0) then
            allocate(A_ref(m, n), B_ref(m_big, nrhs), work_ref(max(1, max(m, n) * 64)))
            A_ref = A_glob
            B_ref = B_glob
            call dgels(transes(ic), m, n, nrhs, A_ref, m, B_ref, m_big, &
                       work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat(B_got, B_ref)
            tol = 64.0_ep * real(max(m, n), ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0,a,i0)') 'trans=', transes(ic), &
                ',m=', m, ',n=', n, ',nrhs=', nrhs
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, work_ref, B_got)
        end if
        deallocate(A_loc, A_glob, B_loc, B_glob, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgels
