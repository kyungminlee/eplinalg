program test_pzgeqpf
    ! Z-mirror of test_pdgeqpf. Same SVD-invariant correctness check.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: zgesvd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix_z, gather_matrix_z
    use target_scalapack, only: target_name, target_eps, target_pzgeqpf
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork, lrwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    complex(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_factored(:,:)
    complex(ep), allocatable :: R(:,:), tau(:), work(:), U(:,:), VT(:,:), A_copy(:,:)
    complex(ep), allocatable :: work_ref(:)
    real(ep),    allocatable :: s_R(:), s_A(:), rwork(:), rwork_ref(:)
    integer,     allocatable :: ipiv(:)
    real(ep) :: err, tol
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1)
    character(len=48) :: label
    integer :: k, j

    call grid_init()
    call report_init('pzgeqpf', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix_z(n, n, mb, nb, A_loc, A_glob, seed = 22801 + 31*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(tau(max(1, locn_a)), ipiv(max(1, locn_a)))
        call target_pzgeqpf(n, n, A_loc, 1, 1, desca, ipiv, tau, &
                            wopt, -1, rwopt, -1, info)
        lwork  = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1)))
        allocate(work(lwork), rwork(lrwork))
        call target_pzgeqpf(n, n, A_loc, 1, 1, desca, ipiv, tau, &
                            work, lwork, rwork, lrwork, info)
        call gather_matrix_z(n, n, mb, nb, A_loc, A_factored)

        if (my_rank == 0) then
            allocate(R(n, n), s_R(n), s_A(n))
            allocate(U(1, 1), VT(1, 1), work_ref(max(1, 64 * n)))
            allocate(rwork_ref(max(1, 5 * n)))
            R = (0.0_ep, 0.0_ep)
            do j = 1, n
                do k = 1, j
                    R(k, j) = A_factored(k, j)
                end do
            end do
            call zgesvd('N', 'N', n, n, R, n, s_R, U, 1, VT, 1, &
                        work_ref, size(work_ref), rwork_ref, info_ref)
            allocate(A_copy(n, n))
            A_copy = A_glob
            call zgesvd('N', 'N', n, n, A_copy, n, s_A, U, 1, VT, 1, &
                        work_ref, size(work_ref), rwork_ref, info_ref)
            err = max_rel_err_vec(s_R, s_A)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(R, s_R, s_A, U, VT, work_ref, rwork_ref, A_copy, A_factored)
        end if
        deallocate(A_loc, A_glob, tau, ipiv, work, rwork)
    end do

    call report_finalize()
    call grid_exit()
end program test_pzgeqpf
