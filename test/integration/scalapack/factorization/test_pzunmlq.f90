program test_pzunmlq
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgelqf, zunmlq
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, &
                                target_pzgelqf, target_pzunmlq
    implicit none

    integer, parameter :: cases = 4
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: transes(*) = ['N', 'C', 'N', 'C']
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, mA, nA, k, mC, nC, info, info_ref, lwork
    integer :: locmA, locnA, lldA, locmC, locnC, lldC
    integer :: desca(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), C_loc(:,:), C_glob(:,:)
    complex(ep), allocatable :: C_got(:,:), A_ref(:,:), C_ref(:,:)
    complex(ep), allocatable :: tau_got(:), tau_ref(:), work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzunmlq', target_name, my_rank)

    do ic = 1, cases
        mA = 32
        if (sides(ic) == 'L') then
            mC = 64; nC = 48
            nA = mC
        else
            mC = 48; nC = 64
            nA = nC
        end if
        k = mA

        call gen_distrib_matrix_z(mA, nA, mb, nb, A_loc, A_glob, seed = 10501 + 17*ic)
        call gen_distrib_matrix_z(mC, nC, mb, nb, C_loc, C_glob, seed = 10511 + 17*ic)

        locmA = numroc_local(mA, mb, my_row, 0, my_nprow)
        locnA = numroc_local(nA, nb, my_col, 0, my_npcol); lldA = max(1, locmA)
        locmC = numroc_local(mC, mb, my_row, 0, my_nprow)
        locnC = numroc_local(nC, nb, my_col, 0, my_npcol); lldC = max(1, locmC)
        call descinit_local(desca, mA, nA, mb, nb, 0, 0, my_context, lldA, info)
        call descinit_local(descc, mC, nC, mb, nb, 0, 0, my_context, lldC, info)

        allocate(tau_got(max(1, locmA)), work(1))
        call target_pzgelqf(mA, nA, A_loc, 1, 1, desca, tau_got, work, -1, info)
        lwork = max(1, int(real(work(1), ep)))
        deallocate(work); allocate(work(lwork))
        call target_pzgelqf(mA, nA, A_loc, 1, 1, desca, tau_got, work, lwork, info)
        deallocate(work)

        allocate(work(1))
        call target_pzunmlq(sides(ic), transes(ic), mC, nC, k, A_loc, 1, 1, desca, &
                            tau_got, C_loc, 1, 1, descc, work, -1, info)
        lwork = max(1, int(real(work(1), ep)))
        deallocate(work); allocate(work(lwork))
        call target_pzunmlq(sides(ic), transes(ic), mC, nC, k, A_loc, 1, 1, desca, &
                            tau_got, C_loc, 1, 1, descc, work, lwork, info)
        call gather_matrix_z(mC, nC, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(A_ref(mA, nA), C_ref(mC, nC), tau_ref(min(mA, nA)), &
                     work_ref(max(1, max(mA, nA, mC, nC) * 64)))
            A_ref = A_glob
            C_ref = C_glob
            call zgelqf(mA, nA, A_ref, mA, tau_ref, work_ref, size(work_ref), info_ref)
            call zunmlq(sides(ic), transes(ic), mC, nC, k, A_ref, mA, tau_ref, &
                        C_ref, mC, work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 64.0_ep * real(max(mA, nA, mC, nC), ep)**2 * target_eps
            write(label, '(a,a,a,a,a,i0,a,i0,a,i0)') 'side=', sides(ic), &
                ',trans=', transes(ic), ',m=', mC, ',n=', nC, ',k=', k
            call report_case(trim(label), err, tol)
            deallocate(A_ref, C_ref, tau_ref, work_ref, C_got)
        end if
        deallocate(A_loc, A_glob, C_loc, C_glob, tau_got, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzunmlq
