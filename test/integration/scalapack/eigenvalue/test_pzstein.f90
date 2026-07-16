program test_pzstein
    ! Complex eigenvectors via inverse iteration on a real symmetric
    ! tridiagonal. Z is complex (zero imaginary part on convergence);
    ! check the residual ||T*v - w*v||_2 / (||T||_inf * ||v||_2) per column.
    use prec_kinds,        only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,   only: dstebz
    use pblas_grid,        only: grid_init, grid_exit, my_rank, my_context, &
                                 my_nprow, my_npcol, my_row, my_col, &
                                 numroc_local, descinit_local
    use pblas_distrib,     only: gather_matrix_z
    use target_scalapack,  only: target_name, target_eps, target_pzstein
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, m, info, nsplit, k
    integer :: locm, locn, lld_z, lwork, liwork
    integer :: descZ(9)
    real(ep), allocatable :: d(:), e(:), W(:), work_st(:), work_dse(:)
    complex(ep), allocatable :: Z_loc(:,:), Z_got(:,:)
    integer,  allocatable :: iblock(:), isplit(:), iwork(:)
    integer,  allocatable :: iwork_st(:), ifail(:), iclustr(:)
    real(ep), allocatable :: gap(:)
    real(ep) :: err, tol, t_norm, max_resid
    character(len=48) :: label

    call grid_init()
    call report_init('pzstein', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)

        allocate(d(n), e(n), W(n), iblock(n), isplit(n))
        block
            integer :: jj, s
            do jj = 1, n
                s = 28401 + 31 * i + jj
                d(jj) = real(mod(s,        1009), ep) / 1009.0_ep + real(2 * jj, ep)
                e(jj) = real(mod(s + 113,  1019), ep) / 1019.0_ep * 0.1_ep
            end do
            e(n) = 0.0_ep
        end block

        allocate(work_dse(4 * n), iwork(3 * n))
        call dstebz('A', 'B', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, &
                    d, e, m, nsplit, W, iblock, isplit, &
                    work_dse, iwork, info)
        deallocate(work_dse, iwork)

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(m, nb, my_col, 0, my_npcol); lld_z = max(1, locm)
        call descinit_local(descZ, n, m, mb, nb, 0, 0, my_context, lld_z, info)
        allocate(Z_loc(lld_z, max(1, locn)))
        Z_loc = (0.0_ep, 0.0_ep)

        allocate(ifail(m), iclustr(2 * my_nprow * my_npcol), gap(my_nprow * my_npcol))
        block
            integer :: p, np00, mq00, iceil_mp
            p = my_nprow * my_npcol
            np00 = numroc_local(n, mb, 0, 0, my_nprow)
            mq00 = numroc_local(m, nb, 0, 0, my_npcol)
            iceil_mp = (m + p - 1) / p
            lwork = max(5 * n, np00 * mq00) + iceil_mp * n + 64
            liwork = 3 * n + p + 1 + 64
        end block
        allocate(work_st(lwork), iwork_st(liwork))

        call target_pzstein(n, d, e, m, W, iblock, isplit, 1.0e-3_ep, &
                            Z_loc, 1, 1, descZ, work_st, lwork, iwork_st, liwork, &
                            ifail, iclustr, gap, info)

        call gather_matrix_z(n, m, mb, nb, Z_loc, Z_got)

        if (my_rank == 0) then
            t_norm = abs(d(1)) + abs(e(1))
            do k = 2, n - 1
                t_norm = max(t_norm, abs(e(k - 1)) + abs(d(k)) + abs(e(k)))
            end do
            t_norm = max(t_norm, abs(e(n - 1)) + abs(d(n)))

            max_resid = 0.0_ep
            block
                integer :: kk, jj
                real(ep) :: v_norm, r
                complex(ep) :: cs
                complex(ep), allocatable :: r_vec(:)
                allocate(r_vec(n))
                do kk = 1, m
                    v_norm = sqrt(sum_sq_z(Z_got(:, kk)))
                    if (v_norm == 0.0_ep) cycle
                    do jj = 1, n
                        cs = d(jj) * Z_got(jj, kk) - W(kk) * Z_got(jj, kk)
                        if (jj > 1) cs = cs + e(jj - 1) * Z_got(jj - 1, kk)
                        if (jj < n) cs = cs + e(jj)     * Z_got(jj + 1, kk)
                        r_vec(jj) = cs
                    end do
                    r = sqrt(sum_sq_z(r_vec))
                    max_resid = max(max_resid, r / (t_norm * v_norm))
                end do
                deallocate(r_vec)
            end block
            err = max_resid
            tol = 64.0_ep * real(n, ep) * target_eps
            write(label, '(a,i0,a,i0)') 'n=', n, ',m=', m
            call report_case(trim(label), err, tol)
            deallocate(Z_got)
        end if
        deallocate(d, e, W, iblock, isplit, Z_loc, work_st, iwork_st, &
                   ifail, iclustr, gap)
    end do

    call report_finalize()
    call grid_exit()

contains
    pure function sum_sq_z(v) result(s)
        complex(ep), intent(in) :: v(:)
        real(ep) :: s
        integer :: jj
        s = 0.0_ep
        do jj = 1, size(v)
            s = s + real(v(jj))**2 + aimag(v(jj))**2
        end do
    end function sum_sq_z
end program test_pzstein
