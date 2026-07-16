program test_pdormtr
    ! Apply Q from sytrd to a matrix C; compare to dsytrd+dormtr.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dsytrd, dormtr
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, &
                                target_pdsytrd, target_pdormtr
    implicit none

    integer, parameter :: cases = 4
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: transes(*) = ['N', 'T', 'N', 'T']
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'L']
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, n, mC, nC, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a, locm_c, locn_c, lld_c
    integer :: desca(9), descc(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), C_loc(:,:), C_glob(:,:)
    real(ep), allocatable :: A_sym(:,:), A_ref(:,:), C_ref(:,:), C_got(:,:)
    real(ep), allocatable :: d(:), e(:), tau(:), work(:)
    real(ep), allocatable :: d_ref(:), e_ref(:), tau_ref(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pdormtr', target_name, my_rank)

    do ic = 1, cases
        n = 64
        if (sides(ic) == 'L') then
            mC = n; nC = 32
        else
            mC = 32; nC = n
        end if

        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 15001 + 17*ic)
        call gen_distrib_matrix(mC, nC, mb, nb, C_loc, C_glob, seed = 15011 + 17*ic)

        allocate(A_sym(n, n))
        A_sym = 0.5_ep * (A_glob + transpose(A_glob))

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_sym(ig, jg)
                end do
            end do
        end if
        locm_c = numroc_local(mC, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(nC, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
        call descinit_local(desca, n,  n,  mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descc, mC, nC, mb, nb, 0, 0, my_context, lld_c, info)

        ! Reduce A to tridiagonal in place.
        allocate(d(max(1, locn_a)), e(max(1, locn_a)), tau(max(1, locn_a)), work(1))
        call target_pdsytrd(uplos(ic), n, A_loc, 1, 1, desca, d, e, tau, work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdsytrd(uplos(ic), n, A_loc, 1, 1, desca, d, e, tau, work, lwork, info)
        deallocate(work)

        ! Apply Q to C.
        allocate(work(1))
        call target_pdormtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, -1, info)
        lwork = max(1, int(work(1)))
        deallocate(work); allocate(work(lwork))
        call target_pdormtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, lwork, info)
        call gather_matrix(mC, nC, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), C_ref(mC, nC), &
                     d_ref(n), e_ref(max(1, n - 1)), tau_ref(max(1, n - 1)), &
                     work_ref(max(1, max(n, mC, nC) * 64)))
            A_ref = A_sym
            C_ref = C_glob
            call dsytrd(uplos(ic), n, A_ref, n, d_ref, e_ref, tau_ref, &
                        work_ref, size(work_ref), info_ref)
            call dormtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                        A_ref, n, tau_ref, C_ref, mC, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat(C_got, C_ref)
            ! pdormtr stitches Householder reflectors with extra
            ! reductions per panel; the L+T path accumulates 4-5 bits
            ! more rounding than U+N. 32*N^3*eps is loose enough not to
            ! flake while still catching real regressions.
            tol = 32.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'side=', sides(ic), &
                ',uplo=', uplos(ic), ',trans=', transes(ic), &
                ',m=', mC, ',n=', nC
            call report_case(trim(label), err, tol)
            deallocate(A_ref, C_ref, d_ref, e_ref, tau_ref, work_ref, C_got)
        end if
        deallocate(A_loc, A_glob, C_loc, C_glob, A_sym, d, e, tau, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdormtr
