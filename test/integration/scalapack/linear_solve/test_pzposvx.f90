program test_pzposvx
    ! Expert Hermitian-positive-definite complex solver. FACT='N' /
    ! EQUED='N': pzposvx Cholesky-factors A internally and solves;
    ! compare gathered X against zposv.
    use prec_kinds,        only: ep
    use compare,           only: max_rel_err_mat_z
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: zposv
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, target_pzposvx
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: nrhs    = 2
    integer, parameter :: mb = 8, nb = 8
    character(len=1), parameter :: uplos(*) = ['U', 'L', 'U']
    integer :: i, n, info, info_ref, lwork, lrwork, k
    integer :: locm_a, locn_a, lld_a, locm_b, locn_b, lld_b
    integer :: desca(9), descaf(9), descb(9), descx(9)
    integer :: owner_r, owner_c, il, jl, ig, jg
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:)
    complex(ep), allocatable :: B_loc(:,:), B_glob(:,:)
    complex(ep), allocatable :: AF_loc(:,:), X_loc(:,:), X_got(:,:)
    complex(ep), allocatable :: A_herm(:,:), A_ref(:,:), B_ref(:,:)
    real(ep), allocatable :: SR(:), SC(:), ferr(:), berr(:)
    complex(ep), allocatable :: work(:)
    real(ep),    allocatable :: rwork(:)
    real(ep) :: err, tol, rcond
    character(len=1) :: equed
    character(len=48) :: label

    call grid_init()
    call report_init('pzposvx', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n,    mb, nb, A_loc, A_glob, seed = 36301 + 31*i)
        call gen_distrib_matrix_z(n, nrhs, mb, nb, B_loc, B_glob, seed = 36311 + 31*i)

        allocate(A_herm(n, n))
        A_herm = 0.5_ep * (A_glob + transpose(conjg(A_glob)))
        do k = 1, n
            A_herm(k, k) = A_herm(k, k) + cmplx(real(2 * n, ep), 0.0_ep, ep)
        end do

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        locm_b = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_b = numroc_local(nrhs, nb, my_col, 0, my_npcol); lld_b = max(1, locm_b)
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
        call descinit_local(desca,  n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descaf, n, n,    mb, nb, 0, 0, my_context, lld_a, info)
        call descinit_local(descb,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)
        call descinit_local(descx,  n, nrhs, mb, nb, 0, 0, my_context, lld_b, info)

        allocate(AF_loc(lld_a, max(1, locn_a)), X_loc(lld_b, max(1, locn_b)))
        AF_loc = (0.0_ep, 0.0_ep); X_loc = (0.0_ep, 0.0_ep)

        allocate(SR(n), SC(n), ferr(nrhs), berr(nrhs))
        SR = 0.0_ep; SC = 0.0_ep
        equed = 'N'

        allocate(work(1), rwork(1))
        call target_pzposvx('N', uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, equed, SR, SC, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, -1, rwork, -1, info)
        lwork  = max(1, int(real(work(1))))
        lrwork = max(1, int(rwork(1)))
        deallocate(work, rwork)
        allocate(work(lwork), rwork(lrwork))
        call target_pzposvx('N', uplos(i), n, nrhs, A_loc, 1, 1, desca, &
                            AF_loc, 1, 1, descaf, equed, SR, SC, &
                            B_loc, 1, 1, descb, X_loc, 1, 1, descx, &
                            rcond, ferr, berr, work, lwork, rwork, lrwork, info)

        call gather_matrix_z(n, nrhs, mb, nb, X_loc, X_got)

        if (my_rank == 0) then
            allocate(A_ref(n, n), B_ref(n, nrhs))
            A_ref = A_herm
            B_ref = B_glob
            call zposv(uplos(i), n, nrhs, A_ref, n, B_ref, n, info_ref)
            err = max_rel_err_mat_z(X_got, B_ref)
            tol = 64.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'uplo=', uplos(i), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, B_ref, X_got)
        end if
        deallocate(A_loc, A_glob, B_loc, B_glob, A_herm, AF_loc, X_loc, &
                   SR, SC, ferr, berr, work, rwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzposvx
