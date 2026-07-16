program test_pdstein
    ! PDSTEIN computes eigenvectors of a symmetric tridiagonal matrix via
    ! inverse iteration, given pre-computed eigenvalues. Inverse iteration
    ! is non-unique up to sign+precision, so we check the residual
    ! ||T*v_i - w_i*v_i|| / (||T||*||v_i||) for each eigenvector instead
    ! of comparing entries directly.
    use prec_kinds,       only: ep
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dstebz
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local, g2l
    use pblas_distrib,    only: gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdstein
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, m, info, nsplit, k
    integer :: locm, locn, lld_z, lwork, liwork
    integer :: descZ(9)
    real(ep), allocatable :: d(:), e(:), W(:)
    real(ep), allocatable :: Z_loc(:,:), Z_got(:,:)
    real(ep), allocatable :: work_st(:), work_dse(:)
    integer,  allocatable :: iblock(:), isplit(:), iwork(:)
    integer,  allocatable :: iwork_st(:), ifail(:), iclustr(:)
    real(ep), allocatable :: gap(:)
    real(ep) :: err, tol, t_norm, max_resid
    character(len=48) :: label

    call grid_init()
    call report_init('pdstein', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)

        allocate(d(n), e(n), W(n), iblock(n), isplit(n))
        block
            integer :: jj, s
            do jj = 1, n
                s = 21401 + 31 * i + jj
                d(jj) = real(mod(s,        1009), ep) / 1009.0_ep + real(2 * jj, ep)
                e(jj) = real(mod(s + 113,  1019), ep) / 1019.0_ep * 0.1_ep
            end do
            e(n) = 0.0_ep
        end block

        ! Replicated bisection on every process to obtain W.
        allocate(work_dse(4 * n), iwork(3 * n))
        call dstebz('A', 'B', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, &
                    d, e, m, nsplit, W, iblock, isplit, &
                    work_dse, iwork, info)
        deallocate(work_dse, iwork)

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(m, nb, my_col, 0, my_npcol); lld_z = max(1, locm)
        call descinit_local(descZ, n, m, mb, nb, 0, 0, my_context, lld_z, info)
        allocate(Z_loc(lld_z, max(1, locn)))
        Z_loc = 0.0_ep

        allocate(ifail(m), iclustr(2 * my_nprow * my_npcol), gap(my_nprow * my_npcol))
        ! LWORK >= max(5*N, NP00*MQ00) + ICEIL(M, P)*N; LIWORK >= 3*N + P + 1
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

        call target_pdstein(n, d, e, m, W, iblock, isplit, 1.0e-3_ep, &
                            Z_loc, 1, 1, descZ, work_st, lwork, iwork_st, liwork, &
                            ifail, iclustr, gap, info)

        call gather_matrix(n, m, mb, nb, Z_loc, Z_got)

        if (my_rank == 0) then
            ! Tridiagonal infinity-norm:  max_i |d_i| + |e_{i-1}| + |e_i|
            t_norm = abs(d(1)) + abs(e(1))
            do k = 2, n - 1
                t_norm = max(t_norm, abs(e(k - 1)) + abs(d(k)) + abs(e(k)))
            end do
            t_norm = max(t_norm, abs(e(n - 1)) + abs(d(n)))

            max_resid = 0.0_ep
            block
                integer :: kk, jj
                real(ep) :: v_norm, r, sum
                real(ep), allocatable :: r_vec(:)
                allocate(r_vec(n))
                do kk = 1, m
                    v_norm = sqrt(sum_sq(Z_got(:, kk)))
                    if (v_norm == 0.0_ep) cycle
                    ! r_j = T_jk * v_k - w*v_j
                    do jj = 1, n
                        sum = d(jj) * Z_got(jj, kk) - W(kk) * Z_got(jj, kk)
                        if (jj > 1) sum = sum + e(jj - 1) * Z_got(jj - 1, kk)
                        if (jj < n) sum = sum + e(jj) * Z_got(jj + 1, kk)
                        r_vec(jj) = sum
                    end do
                    r = sqrt(sum_sq(r_vec))
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
    pure function sum_sq(v) result(s)
        real(ep), intent(in) :: v(:)
        real(ep) :: s
        integer :: jj
        s = 0.0_ep
        do jj = 1, size(v)
            s = s + v(jj) * v(jj)
        end do
    end function sum_sq
end program test_pdstein
