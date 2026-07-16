program test_pdstedc
    ! Divide-and-conquer eigensolver for symmetric tridiagonal matrices.
    ! D becomes eigenvalues; Q (initialized to identity via COMPZ='I')
    ! becomes eigenvectors. Compare eigenvalues against quad dstebz, and
    ! verify residual ||T*v_i - lambda_i*v_i|| / (||T|| * ||v_i||).
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dstebz
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdstedc
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, m, nsplit
    integer :: locm, locn, lld_q, lwork, liwork
    integer :: descq(9)
    real(ep), allocatable :: d(:), e(:), d_in(:), e_in(:), W_ref(:)
    real(ep), allocatable :: Q_loc(:,:), Q_got(:,:)
    real(ep), allocatable :: work(:), work_dse(:)
    integer,  allocatable :: iwork(:), iblock(:), isplit(:), iwork_dse(:)
    real(ep) :: err, tol, t_norm, max_resid
    character(len=48) :: label

    call grid_init()
    call report_init('pdstedc', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        allocate(d(n), e(n), d_in(n), e_in(n))
        block
            integer :: jj, s
            do jj = 1, n
                s = 23301 + 31 * i + jj
                d(jj) = real(mod(s,        1009), ep) / 1009.0_ep + real(2 * jj, ep)
                e(jj) = real(mod(s + 113,  1019), ep) / 1019.0_ep * 0.1_ep
            end do
            e(n) = 0.0_ep
        end block
        d_in = d; e_in = e

        locm = numroc_local(n, mb, my_row, 0, my_nprow)
        locn = numroc_local(n, nb, my_col, 0, my_npcol); lld_q = max(1, locm)
        call descinit_local(descq, n, n, mb, nb, 0, 0, my_context, lld_q, info)
        allocate(Q_loc(lld_q, max(1, locn)))
        Q_loc = 0.0_ep

        ! Workspace query
        allocate(work(1), iwork(1))
        call target_pdstedc('I', n, d, e, Q_loc, 1, 1, descq, &
                            work, -1, iwork, -1, info)
        lwork  = max(1, int(work(1)))
        liwork = max(1, iwork(1))
        deallocate(work, iwork)
        allocate(work(lwork), iwork(liwork))
        ! Reset d/e since the query may have touched them
        d = d_in; e = e_in
        call target_pdstedc('I', n, d, e, Q_loc, 1, 1, descq, &
                            work, lwork, iwork, liwork, info)
        call gather_matrix(n, n, mb, nb, Q_loc, Q_got)

        if (my_rank == 0) then
            ! Reference eigenvalues via quad dstebz
            allocate(W_ref(n), iblock(n), isplit(n))
            allocate(work_dse(4 * n), iwork_dse(3 * n))
            call dstebz('A', 'E', n, 0.0_ep, 0.0_ep, 0, 0, 0.0_ep, &
                        d_in, e_in, m, nsplit, W_ref, iblock, isplit, &
                        work_dse, iwork_dse, info)
            ! pdstedc returns eigenvalues in d (ascending order).
            err = max_rel_err_vec(d(1:n), W_ref(1:n))
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0,a)') 'n=', n, ',out=W'
            call report_case(trim(label), err, tol)

            ! Eigenvector residual.
            t_norm = abs(d_in(1)) + abs(e_in(1))
            block
                integer :: kk
                do kk = 2, n - 1
                    t_norm = max(t_norm, abs(e_in(kk - 1)) + abs(d_in(kk)) + abs(e_in(kk)))
                end do
                t_norm = max(t_norm, abs(e_in(n - 1)) + abs(d_in(n)))
            end block

            max_resid = 0.0_ep
            block
                integer :: kk, jj
                real(ep) :: v_norm, r, sum
                real(ep), allocatable :: r_vec(:)
                allocate(r_vec(n))
                do kk = 1, n
                    v_norm = sqrt(sum_sq(Q_got(:, kk)))
                    if (v_norm == 0.0_ep) cycle
                    do jj = 1, n
                        sum = d_in(jj) * Q_got(jj, kk) - d(kk) * Q_got(jj, kk)
                        if (jj > 1) sum = sum + e_in(jj - 1) * Q_got(jj - 1, kk)
                        if (jj < n) sum = sum + e_in(jj) * Q_got(jj + 1, kk)
                        r_vec(jj) = sum
                    end do
                    r = sqrt(sum_sq(r_vec))
                    max_resid = max(max_resid, r / (t_norm * v_norm))
                end do
                deallocate(r_vec)
            end block
            err = max_resid
            write(label, '(a,i0,a)') 'n=', n, ',out=residual'
            call report_case(trim(label), err, tol)
            deallocate(W_ref, iblock, isplit, work_dse, iwork_dse, Q_got)
        end if
        deallocate(d, e, d_in, e_in, Q_loc, work, iwork)
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
end program test_pdstedc
