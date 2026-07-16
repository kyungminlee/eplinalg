program test_pzunmhr
    ! Apply Q from gehrd to a matrix C; compare to zgehrd+zunmhr.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgehrd, zunmhr
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, &
                                target_pzgehrd, target_pzunmhr
    implicit none

    integer, parameter :: cases = 4
    character(len=1), parameter :: sides(*)   = ['L', 'L', 'R', 'R']
    character(len=1), parameter :: transes(*) = ['N', 'C', 'N', 'C']
    integer, parameter :: mb = 8, nb = 8
    integer :: ic, n, mC, nC, ilo, ihi, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a, locm_c, locn_c, lld_c
    integer :: desca(9), descc(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), C_loc(:,:), C_glob(:,:)
    complex(ep), allocatable :: A_ref(:,:), C_ref(:,:), C_got(:,:)
    complex(ep), allocatable :: tau(:), tau_ref(:), work(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzunmhr', target_name, my_rank)

    do ic = 1, cases
        n = 64; ilo = 1; ihi = n
        if (sides(ic) == 'L') then
            mC = n; nC = 32
        else
            mC = 32; nC = n
        end if

        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 15201 + 17*ic)
        call gen_distrib_matrix_z(mC, nC, mb, nb, C_loc, C_glob, seed = 15211 + 17*ic)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_c = numroc_local(mC, mb, my_row, 0, my_nprow)
        locn_c = numroc_local(nC, nb, my_col, 0, my_npcol); lld_c = max(1, locm_c)
        call descinit_local(desca, n,  n,  mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descc, mC, nC, mb, nb, 0, 0, my_context, lld_c, info)

        allocate(tau(max(1, locn_a)), work(1))
        call target_pzgehrd(n, ilo, ihi, A_loc, 1, 1, desca, tau, work, -1, info)
        lwork = max(1, int(real(work(1), ep)))
        deallocate(work); allocate(work(lwork))
        call target_pzgehrd(n, ilo, ihi, A_loc, 1, 1, desca, tau, work, lwork, info)
        deallocate(work)

        allocate(work(1))
        call target_pzunmhr(sides(ic), transes(ic), mC, nC, ilo, ihi, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, -1, info)
        lwork = max(1, int(real(work(1), ep)))
        deallocate(work); allocate(work(lwork))
        call target_pzunmhr(sides(ic), transes(ic), mC, nC, ilo, ihi, &
                            A_loc, 1, 1, desca, tau, C_loc, 1, 1, descc, &
                            work, lwork, info)
        call gather_matrix_z(mC, nC, mb, nb, C_loc, C_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), C_ref(mC, nC), &
                     tau_ref(max(1, n - 1)), &
                     work_ref(max(1, max(n, mC, nC) * 64)))
            A_ref = A_glob
            C_ref = C_glob
            call zgehrd(n, ilo, ihi, A_ref, n, tau_ref, work_ref, &
                        size(work_ref), info_ref)
            call zunmhr(sides(ic), transes(ic), mC, nC, ilo, ihi, &
                        A_ref, n, tau_ref, C_ref, mC, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(C_got, C_ref)
            tol = 32.0_ep * real(n, ep)**3 * target_eps
            write(label, '(a,a,a,a,a,i0,a,i0)') 'side=', sides(ic), &
                ',trans=', transes(ic), ',m=', mC, ',n=', nC
            call report_case(trim(label), err, tol)
            deallocate(A_ref, C_ref, tau_ref, work_ref, C_got)
        end if
        deallocate(A_loc, A_glob, C_loc, C_glob, tau, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzunmhr
