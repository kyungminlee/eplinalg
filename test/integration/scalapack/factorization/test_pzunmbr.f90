program test_pzunmbr
    ! Apply Q (or P^H) from gebrd to a matrix C; compare to zgebrd+zunmbr.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: zgebrd, zunmbr
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local
    use pblas_distrib,     only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pzgebrd, target_pzunmbr
    implicit none

    ! See test_pdormbr.f90 — only VECT='Q' branches are exercised pending
    ! a descA/descC alignment study for VECT='P'.
    integer, parameter :: cases = 2
    character(len=1), parameter :: vects(*)   = ['Q', 'Q']
    character(len=1), parameter :: sides(*)   = ['L', 'R']
    character(len=1), parameter :: transes(*) = ['N', 'C']
    integer, parameter :: m_a = 64, n_a = 48
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, mC, nC, k_arg, info, info_ref, lwork, mn
    integer :: locm_a, locn_a, lld_a, locm_c, locn_c, lld_c
    integer :: desca(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), C_loc(:,:), C_glob(:,:)
    complex(ep), allocatable :: A_ref(:,:), C_ref(:,:), C_got(:,:)
    real(ep),    allocatable :: d(:), e(:)
    complex(ep), allocatable :: tauq(:), taup(:), work(:)
    real(ep),    allocatable :: d_ref(:), e_ref(:)
    complex(ep), allocatable :: tauq_ref(:), taup_ref(:), work_ref(:)
    real(ep) :: err, tol
    character(len=64) :: label

    call grid_init()
    call report_init('pzunmbr', target_name, my_rank)

    mn = min(m_a, n_a)

    call gen_distrib_matrix_z(m_a, n_a, mb, nb, A_loc, A_glob, seed = 24501)

    locm_a = numroc_local(m_a, mb, my_row, 0, my_nprow)
    locn_a = numroc_local(n_a, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
    call descinit_local(desca, m_a, n_a, mb, nb, 0, 0, my_context, lld_a, info)

    allocate(d(max(1,locn_a)), e(max(1,locm_a)), &
             tauq(max(1,locn_a)), taup(max(1,locm_a)), work(1))
    call target_pzgebrd(m_a, n_a, A_loc, 1, 1, desca, d, e, tauq, taup, &
                        work, -1, info)
    lwork = max(1, int(real(work(1))))
    deallocate(work); allocate(work(lwork))
    call target_pzgebrd(m_a, n_a, A_loc, 1, 1, desca, d, e, tauq, taup, &
                        work, lwork, info)
    deallocate(work)

    if (my_rank == 0) then
        allocate(A_ref(m_a, n_a), d_ref(mn), e_ref(mn), &
                 tauq_ref(mn), taup_ref(mn), &
                 work_ref(max(1, (m_a + n_a) * 64)))
        A_ref = A_glob
        call zgebrd(m_a, n_a, A_ref, m_a, d_ref, e_ref, tauq_ref, taup_ref, &
                    work_ref, size(work_ref), info_ref)
    end if

    do ic = 1, cases
        if (sides(ic) == 'L') then; mC = m_a; nC = 24
        else;                       mC = 24;  nC = m_a; end if
        k_arg = mn

        call gen_distrib_matrix_z(mC, nC, mb, nb, C_loc, C_glob, seed = 24601 + 13*ic)

        locm_c = numroc_local(mC, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(nC, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
        call descinit_local(descc, mC, nC, mb, nb, 0, 0, my_context, lld_c, info)

        allocate(work(1))
        call target_pzunmbr(vects(ic), sides(ic), transes(ic), &
                            mC, nC, k_arg, A_loc, 1, 1, desca, tauq, &
                            C_loc, 1, 1, descc, work, -1, info)
        lwork = max(1, int(real(work(1))))
        deallocate(work); allocate(work(lwork))
        call target_pzunmbr(vects(ic), sides(ic), transes(ic), &
                            mC, nC, k_arg, A_loc, 1, 1, desca, tauq, &
                            C_loc, 1, 1, descc, work, lwork, info)
        call gather_matrix_z(mC, nC, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(C_ref(mC, nC))
            C_ref = C_glob
            call zunmbr(vects(ic), sides(ic), transes(ic), &
                        mC, nC, k_arg, A_ref, m_a, tauq_ref, &
                        C_ref, mC, work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 32.0_ep * real(max(m_a, n_a), ep)**3 * target_eps
            write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'vect=', vects(ic), &
                ',side=', sides(ic), ',trans=', transes(ic), &
                ',m=', mC, ',n=', nC
            call report_case(trim(label), err, tol)
            deallocate(C_ref, C_got)
        end if
        deallocate(C_loc, C_glob, work)
    end do

    if (my_rank == 0) deallocate(A_ref, d_ref, e_ref, tauq_ref, taup_ref, work_ref)
    deallocate(A_loc, A_glob, d, e, tauq, taup)

    call report_finalize()
    call grid_exit()
end program test_pzunmbr
