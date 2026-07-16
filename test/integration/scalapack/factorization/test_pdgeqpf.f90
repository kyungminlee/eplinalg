program test_pdgeqpf
    ! QR with column pivoting: A*P = Q*R, where P is the column
    ! permutation given by IPIV. Pivot ordering can diverge between
    ! parallel and serial paths, so we don't compare IPIV directly.
    ! Instead, use the invariant that SVD(A) = SVD(R) (up to ordering),
    ! since A*P = Q*R and Q, P are orthogonal/permutation.
    use prec_kinds,       only: ep
    use compare,          only: max_rel_err_vec
    use pblas_prec_report, only: report_init, report_case, report_finalize
    use ref_quad_lapack,  only: dgesvd
    use pblas_grid,       only: grid_init, grid_exit, my_rank, my_context, &
                                my_nprow, my_npcol, my_row, my_col, &
                                numroc_local, descinit_local
    use pblas_distrib,    only: gen_distrib_matrix, gather_matrix
    use target_scalapack, only: target_name, target_eps, target_pdgeqpf
    implicit none

    integer, parameter :: ns(*) = [32, 64, 96]
    integer, parameter :: mb = 8, nb = 8
    integer :: i, n, info, info_ref, lwork
    integer :: locm_a, locn_a, lld_a
    integer :: desca(9)
    real(ep), allocatable :: A_loc(:,:), A_glob(:,:), A_factored(:,:)
    real(ep), allocatable :: R(:,:), tau(:), work(:)
    real(ep), allocatable :: s_R(:), s_A(:), work_ref(:), U(:,:), VT(:,:)
    real(ep), allocatable :: A_copy(:,:)
    integer,  allocatable :: ipiv(:)
    real(ep) :: err, tol
    real(ep) :: wopt(1)
    character(len=48) :: label
    integer :: k, j

    call grid_init()
    call report_init('pdgeqpf', target_name, my_rank)

    do i = 1, size(ns)
        n = ns(i)
        call gen_distrib_matrix(n, n, mb, nb, A_loc, A_glob, seed = 22701 + 31*i)

        locm_a = numroc_local(n, mb, my_row, 0, my_nprow)
        locn_a = numroc_local(n, nb, my_col, 0, my_npcol); lld_a = max(1, locm_a)
        call descinit_local(desca, n, n, mb, nb, 0, 0, my_context, lld_a, info)

        allocate(tau(max(1, locn_a)), ipiv(max(1, locn_a)))
        call target_pdgeqpf(n, n, A_loc, 1, 1, desca, ipiv, tau, &
                            wopt, -1, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call target_pdgeqpf(n, n, A_loc, 1, 1, desca, ipiv, tau, &
                            work, lwork, info)
        call gather_matrix(n, n, mb, nb, A_loc, A_factored)

        if (my_rank == 0) then
            allocate(R(n, n), s_R(n), s_A(n))
            allocate(U(1, 1), VT(1, 1), work_ref(max(1, 64 * n)))
            R = 0.0_ep
            do j = 1, n
                do k = 1, j
                    R(k, j) = A_factored(k, j)
                end do
            end do
            ! SVD of R
            call dgesvd('N', 'N', n, n, R, n, s_R, U, 1, VT, 1, &
                        work_ref, size(work_ref), info_ref)
            ! SVD of original A
            allocate(A_copy(n, n))
            A_copy = A_glob
            call dgesvd('N', 'N', n, n, A_copy, n, s_A, U, 1, VT, 1, &
                        work_ref, size(work_ref), info_ref)
            err = max_rel_err_vec(s_R, s_A)
            tol = 32.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,i0)') 'n=', n
            call report_case(trim(label), err, tol)
            deallocate(R, s_R, s_A, U, VT, work_ref, A_copy, A_factored)
        end if
        deallocate(A_loc, A_glob, tau, ipiv, work)
    end do

    call report_finalize()
    call grid_exit()
end program test_pdgeqpf
