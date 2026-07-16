program test_pdgesvd
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgesvd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgesvd
    implicit none

    integer, parameter :: ms(*) = [48, 80, 96]
    integer, parameter :: ns(*) = [32, 48, 64]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, m, n, info, info_ref, lwork, minmn
    integer :: locm_a, locn_a, lld_a, locm_u, locn_vt, lld_u, lld_vt
    integer :: desca(9), descu(9), descvt(9)
    real(ep), allocatable :: A_loc(:,:), U_loc(:,:), VT_loc(:,:)
    real(ep), allocatable :: A_glob(:,:), A_ref(:,:)
    real(ep), allocatable :: s_got(:), s_ref(:), work(:), work_ref(:)
    real(ep), allocatable :: Uref(:,:), VTref(:,:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdgesvd', target_name, my_rank)

    do i = 1, size(ms)
        m = ms(i); n = ns(i); minmn = min(m, n)
        call gen_distrib_matrix(m, n, mb, nb, A_loc, A_glob, seed = 9801 + 31*i)

        locm_a  = numroc_local(m, mb, my_row, 0, my_nprow)
        locn_a  = numroc_local(n, nb, my_col, 0, my_npcol); lld_a  = max(1, locm_a)
        locm_u  = numroc_local(m,     mb, my_row, 0, my_nprow); lld_u  = max(1, locm_u)
        locn_vt = numroc_local(n,     nb, my_col, 0, my_npcol)
        ! VT is a minmn × n descriptor; its LLD is row-block-cyclic
        ! over minmn, not n. The two coincide when m ≥ n (current test
        ! shapes), but the wider m < n case would silently corrupt the
        ! VT layout if this used n.
        lld_vt  = max(1, numroc_local(minmn, mb, my_row, 0, my_nprow))
        call descinit_local(desca,  m, n, mb, nb, 0, 0, my_context, lld_a,  info)
        call descinit_local(descu,  m, minmn, mb, nb, 0, 0, my_context, lld_u,  info)
        call descinit_local(descvt, minmn, n, mb, nb, 0, 0, my_context, lld_vt, info)

        allocate(U_loc(lld_u, max(1, numroc_local(minmn, nb, my_col, 0, my_npcol))))
        allocate(VT_loc(lld_vt, max(1, locn_vt)))
        U_loc  = 0.0_ep
        VT_loc = 0.0_ep

        ! Workspace query — same approach as pdsyev: PDGESVD's exact
        ! requirement depends on minmn, block-cyclic locals, and the
        ! BLACS grid shape, so let the routine itself report it.
        allocate(s_got(minmn), work(1))
        call target_pdgesvd('V', 'V', m, n, A_loc, 1, 1, desca, s_got, &
                            U_loc, 1, 1, descu, VT_loc, 1, 1, descvt, &
                            work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work)
        allocate(work(lwork))
        call target_pdgesvd('V', 'V', m, n, A_loc, 1, 1, desca, s_got, &
                            U_loc, 1, 1, descu, VT_loc, 1, 1, descvt, &
                            work, lwork, info)

        if (my_rank == 0) then
            allocate(A_ref(m, n), s_ref(minmn))
            allocate(Uref(m, minmn), VTref(minmn, n))
            allocate(work_ref(max(1, 64 * max(m, n))))
            A_ref = A_glob
            call dgesvd('S', 'S', m, n, A_ref, m, s_ref, Uref, m, &
                        VTref, minmn, work_ref, size(work_ref), info_ref)
            err = max_rel_err_vec(s_got, s_ref)
            tol = 32.0_ep * real(max(m, n), ep)**2 * target_eps
            write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, s_ref, Uref, VTref, work_ref)
        end if
        deallocate(A_loc, U_loc, VT_loc, A_glob, s_got, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgesvd
