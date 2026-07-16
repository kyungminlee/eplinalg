program test_pzunmtr
    ! Apply Q from hetrd to a matrix C; compare to zhetrd+zunmtr.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: zhetrd, zunmtr
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, &
                                 target_pzhetrd, target_pzunmtr
    implicit none

    integer, parameter :: cases = 4
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: transes(*) = ['N', 'C', 'N', 'C']
    character(len=1), parameter :: uplos(*)   = ['U', 'L', 'U', 'L']
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, n, mC, nC, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a, locm_c, locn_c, lld_c
    integer :: desca(9), descc(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), C_loc(:,:), C_glob(:,:)
    complex(ep), allocatable :: A_herm(:,:), A_ref(:,:), C_ref(:,:), C_got(:,:)
    real(ep),    allocatable :: d(:), e(:)
    complex(ep), allocatable :: tau(:), work(:)
    real(ep),    allocatable :: d_ref(:), e_ref(:)
    complex(ep), allocatable :: tau_ref(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzunmtr', target_name, my_rank)

    do ic = 1, cases
        n = 64
        if (sides(ic) == 'L') then
            mC = n; nC = 32
        else
            mC = 32; nC = n
        end if

        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 25001 + 17*ic)
        call gen_distrib_matrix_z(mC, nC, mb, nb, C_loc, C_glob, seed = 25011 + 17*ic)

        allocate(A_herm(n, n))
        A_herm = 0.5_ep * (A_glob + conjg(transpose(A_glob)))

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        if (locm_a > 0 .and. locn_a > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) A_loc(il, jl) = A_herm(ig, jg)
                end do
            end do
        end if
        locm_c = numroc_local(mC, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(nC, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
        call descinit_local(desca, n,  n,  mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descc, mC, nC, mb, nb, 0, 0, my_context, lld_c, info)

        allocate(d(max(1, locn_a)), e(max(1, locn_a)), &
                 tau(max(1, locn_a)), work(1))
        call target_pzhetrd(uplos(ic), n, A_loc, 1, 1, desca, d, e, tau, &
                            work, -1, info)
        lwork = max(1, int(real(work(1))))
        deallocate(work); allocate(work(lwork))
        call target_pzhetrd(uplos(ic), n, A_loc, 1, 1, desca, d, e, tau, &
                            work, lwork, info)
        deallocate(work)

        allocate(work(1))
        call target_pzunmtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, -1, info)
        lwork = max(1, int(real(work(1))))
        deallocate(work); allocate(work(lwork))
        call target_pzunmtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, lwork, info)
        call gather_matrix_z(mC, nC, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), C_ref(mC, nC), &
                     d_ref(n), e_ref(max(1, n - 1)), tau_ref(max(1, n - 1)), &
                     work_ref(max(1, max(n, mC, nC) * 64)))
            A_ref = A_herm
            C_ref = C_glob
            call zhetrd(uplos(ic), n, A_ref, n, d_ref, e_ref, tau_ref, &
                        work_ref, size(work_ref), info_ref)
            call zunmtr(sides(ic), uplos(ic), transes(ic), mC, nC, &
                        A_ref, n, tau_ref, C_ref, mC, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 32.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,a,a,a,a,i0,a,i0)') 'side=', sides(ic), &
                ',uplo=', uplos(ic), ',trans=', transes(ic), &
                ',m=', mC, ',n=', nC
            call report_case(trim(label), err, tol)
            deallocate(A_ref, C_ref, d_ref, e_ref, tau_ref, work_ref, C_got)
        end if
        deallocate(A_loc, A_glob, C_loc, C_glob, A_herm, d, e, tau, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzunmtr
