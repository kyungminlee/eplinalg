program test_pztrevc
    ! Right-eigenvectors of a complex upper-triangular matrix T.
    ! For an upper-triangular T, the eigenvalues are the diagonal
    ! entries; eigenvectors are well-defined modulo scale + phase.
    ! Verify by checking the residual ||T*v_k - lambda_k * v_k|| is
    ! small for each computed VR(:, k).
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local, g2l
    use pblas_distrib,     only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, target_pztrevc
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, lwork, lrwork, m, k, j
    integer :: locm, locn, lld
    integer :: descT(9), descVL(9), descVR(9)
    integer :: ig, jg, owner_r, owner_c, il, jl
    complex(ep), allocatable :: T_loc(:,:), T_glob(:,:)
    complex(ep), allocatable :: VL_loc(:,:), VR_loc(:,:), VR_got(:,:)
    complex(ep), allocatable :: work(:)
    real(ep),    allocatable :: rwork(:)
    logical,     allocatable :: select_v(:)
    real(ep) :: err, tol, t_norm, max_resid, v_norm, r
    complex(ep) :: lambda, cs
    complex(ep), allocatable :: r_vec(:)
    character(len=48) :: label

    call grid_init()
    call report_init('pztrevc', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, T_loc, T_glob, seed = 32101 + 31*i)

        do jg = 1, n
            do ig = jg + 1, n
                T_glob(ig, jg) = (0.0_ep, 0.0_ep)
            end do
            ! Bias diagonal to keep eigenvalues well-separated and large.
            T_glob(jg, jg) = T_glob(jg, jg) + cmplx(real(2 * jg, ep), 0.0_ep, ep)
        end do
        if (size(T_loc, 1) > 0 .and. size(T_loc, 2) > 0) then
            do jg = 1, n
                call g2l(jg, nb, my_npcol, owner_c, jl)
                if (owner_c /= my_col) cycle
                do ig = 1, n
                    call g2l(ig, mb, my_nprow, owner_r, il)
                    if (owner_r == my_row) T_loc(il, jl) = T_glob(ig, jg)
                end do
            end do
        end if

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol); lld = max(1, locm)
        call descinit_local(descT,  n, n, mb, nb, 0, 0, my_context, lld, info)
        call descinit_local(descVL, n, n, mb, nb, 0, 0, my_context, lld, info)
        call descinit_local(descVR, n, n, mb, nb, 0, 0, my_context, lld, info)

        allocate(VL_loc(lld, max(1, locn)), VR_loc(lld, max(1, locn)))
        VL_loc = (0.0_ep, 0.0_ep); VR_loc = (0.0_ep, 0.0_ep)
        allocate(select_v(n)); select_v = .true.

        lwork  = 2 * n + 64; lrwork = n + 64
        allocate(work(lwork), rwork(lrwork))
        call target_pztrevc('R', 'A', select_v, n, T_loc, descT, &
                            VL_loc, descVL, VR_loc, descVR, n, m, &
                            work, rwork, info)

        call gather_matrix_z(n, n, mb, nb, VR_loc, VR_got)

        if (my_rank == 0) then
            ! Frobenius norm of T (use global T_glob, which is replicated).
            t_norm = 0.0_ep
            do j = 1, n
                do k = 1, j
                    t_norm = t_norm + real(T_glob(k, j))**2 + aimag(T_glob(k, j))**2
                end do
            end do
            t_norm = sqrt(t_norm)

            allocate(r_vec(n))
            max_resid = 0.0_ep
            do k = 1, m
                lambda = T_glob(k, k)
                v_norm = 0.0_ep
                do j = 1, n
                    v_norm = v_norm + real(VR_got(j, k))**2 + aimag(VR_got(j, k))**2
                end do
                v_norm = sqrt(v_norm)
                if (v_norm == 0.0_ep) cycle
                ! r = T*v - lambda*v
                do j = 1, n
                    cs = (0.0_ep, 0.0_ep)
                    do ig = 1, j
                        cs = cs + T_glob(ig, j) * VR_got(ig, k)  ! upper triangle: row contribution
                    end do
                    ! correction: above accumulates T(:,j)*v(j) wrong; rewrite
                end do
                ! Correct: r_i = sum_p T(i,p) v(p,k) - lambda * v(i,k), p>=i
                do ig = 1, n
                    cs = (0.0_ep, 0.0_ep)
                    do j = ig, n
                        cs = cs + T_glob(ig, j) * VR_got(j, k)
                    end do
                    cs = cs - lambda * VR_got(ig, k)
                    r_vec(ig) = cs
                end do
                r = 0.0_ep
                do ig = 1, n
                    r = r + real(r_vec(ig))**2 + aimag(r_vec(ig))**2
                end do
                r = sqrt(r)
                max_resid = max(max_resid, r / (t_norm * v_norm))
            end do
            deallocate(r_vec)
            err = max_resid
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',m=', m
            call report_case(trim(label), err, tol)
            deallocate(VR_got)
        end if
        deallocate(T_loc, T_glob, VL_loc, VR_loc, select_v, work, rwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pztrevc
