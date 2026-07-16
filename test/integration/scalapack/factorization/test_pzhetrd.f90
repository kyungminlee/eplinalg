program test_pzhetrd
    ! Hermitian tridiag reduction A -> Q*T*Q^H. Compare gathered A
    ! (encodes T + Householder reflectors) against zhetrd.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zhetrd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzhetrd
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_got(:,:)
    complex(ep), allocatable :: A_herm(:,:), A_ref(:,:)
    real(ep), allocatable :: d(:), e(:)
    complex(ep), allocatable :: tau(:), work(:)
    real(ep), allocatable :: d_ref(:), e_ref(:)
    complex(ep), allocatable :: tau_ref(:), work_ref(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call grid_init()
    call report_init('pzhetrd', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 14801 + 31*i)

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
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(d(max(1, locn_a)), e(max(1, locn_a)), tau(max(1, locn_a)), work(1))
        call target_pzhetrd(uplos(i), n, A_loc, 1, 1, desca, d, e, tau, &
                            work, -1, info)
        lwork = max(1, int(real(work(1), ep)))
        deallocate(work); allocate(work(lwork))
        call target_pzhetrd(uplos(i), n, A_loc, 1, 1, desca, d, e, tau, &
                            work, lwork, info)
        call gather_matrix_z(n, n, mb, nb, A_loc, A_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), d_ref(n), e_ref(max(1, n - 1)), &
                     tau_ref(max(1, n - 1)), work_ref(max(1, n * 64)))
            A_ref = A_herm
            call zhetrd(uplos(i), n, A_ref, n, d_ref, e_ref, tau_ref, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_mat_z(A_got, A_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, d_ref, e_ref, tau_ref, work_ref, A_got)
        end if
        deallocate(A_loc, A_glob, A_herm, d, e, tau, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzhetrd
