program test_pzheevr
    ! Hermitian eigensolver via Relatively Robust Representations.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zheev
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzheevr
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, lrwork, liwork, m, nz, j
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9), descz(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: A_loc(:,:), Z_loc(:,:), Z_glob(:,:)
    complex(ep), allocatable :: A_glob(:,:), A_herm(:,:), A_ref(:,:)
    complex(ep), allocatable :: work(:), work_ref(:)
    real(ep),    allocatable :: w_got(:), w_ref(:), rwork(:), rwork_ref(:)
    integer,     allocatable :: iwork(:)
    complex(ep), allocatable :: AZ(:,:), ZD(:,:)
    real(ep) :: err, tol, anorm, rnorm
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1)
    integer  :: iwopt(1)

    call grid_init()
    call report_init('pzheevr', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 21701 + 31*i)
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
        call descinit_local(descz, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(Z_loc(max(1, locm_a), max(1, locn_a)))
        Z_loc = (0.0_ep, 0.0_ep)
        allocate(w_got(n))

        call target_pzheevr('V', 'A', 'U', n, A_loc, 1, 1, desca, &
                            0.0_ep, 0.0_ep, 0, 0, m, nz, w_got, &
                            Z_loc, 1, 1, descz, wopt, -1, rwopt, -1, &
                            iwopt, -1, info)
        lwork  = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        call target_pzheevr('V', 'A', 'U', n, A_loc, 1, 1, desca, &
                            0.0_ep, 0.0_ep, 0, 0, m, nz, w_got, &
                            Z_loc, 1, 1, descz, work, lwork, rwork, lrwork, &
                            iwork, liwork, info)
        call gather_matrix_z(n, n, mb, nb, Z_loc, Z_glob)

        if (my_rank == 0) then
            allocate(A_ref(n, n), w_ref(n), work_ref(max(1, n * 64)), &
                     rwork_ref(max(1, 3 * n - 2)))
            A_ref = A_herm
            call zheev('V', 'U', n, A_ref, n, w_ref, &
                       work_ref, size(work_ref), rwork_ref, info_ref)
            err = max_rel_err_vec(w_got, w_ref)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            block
                character(len=48) :: label
                write(label, '(a,i0,a)') 'n=', n, ',out=W'
                call report_case(trim(label), err, tol)
            end block
            allocate(AZ(n, n), ZD(n, n))
            AZ = matmul(A_herm, Z_glob)
            do j = 1, n
                ZD(:, j) = Z_glob(:, j) * w_got(j)
            end do
            anorm = maxval(abs(A_herm))
            rnorm = maxval(abs(AZ - ZD))
            err = rnorm / max(anorm, tiny(1.0_ep))
            block
                character(len=48) :: label
                write(label, '(a,i0,a)') 'n=', n, ',out=residual'
                call report_case(trim(label), err, tol)
            end block
            deallocate(AZ, ZD, A_ref, w_ref, work_ref, rwork_ref)
        end if
        if (allocated(Z_glob)) deallocate(Z_glob)
        deallocate(A_loc, Z_loc, A_glob, A_herm, w_got, work, rwork, iwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzheevr
