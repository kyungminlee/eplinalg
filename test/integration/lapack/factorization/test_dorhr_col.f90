! dorhr_col: Householder reconstruction from Q (column variant).
program test_dorhr_col
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dorhr_col
    use ref_quad_lapack, only: dlatsqr, dorgtsqr, dorhr_col
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [4,  8]
    integer, parameter :: mbs(*) = [8, 16]
    integer, parameter :: nbs(*) = [2, 4]
    integer :: i, m, n, mb, nb, info, lwork, ntcols
    real(ep), allocatable :: A0(:,:), Q(:,:), A_ref(:,:), A_got(:,:)
    real(ep), allocatable :: T_arr(:,:), T_ref(:,:), T_got(:,:)
    real(ep), allocatable :: D_ref(:), D_got(:), work(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorhr_col', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); mb = mbs(i); nb = nbs(i)
        call gen_matrix_quad(m, n, A0, seed = 23401 + 89 * i)
        ntcols = n * max(1, (m - n + (mb - n) - 1) / max(1, mb - n))
        allocate(Q(m, n), T_arr(nb, ntcols))
        Q = A0
        call dlatsqr(m, n, mb, nb, Q, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dlatsqr(m, n, mb, nb, Q, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        call dorgtsqr(m, n, mb, nb, Q, m, T_arr, nb, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorgtsqr(m, n, mb, nb, Q, m, T_arr, nb, work, lwork, info)
        deallocate(work)
        ! Now Q is an orthonormal m x n matrix. Reconstruct via dorhr_col.
        allocate(A_ref(m, n), A_got(m, n), T_ref(nb, n), T_got(nb, n), D_ref(n), D_got(n))
        A_ref = Q; A_got = Q
        T_ref = 0.0_ep; T_got = 0.0_ep
        call dorhr_col(m, n, nb, A_ref, m, T_ref, nb, D_ref, info)
        call target_dorhr_col(m, n, nb, A_got, m, T_got, nb, D_got, info)
        err = max(max_rel_err_mat(A_got, A_ref), max_rel_err_mat(T_got, T_ref), &
                  max_rel_err_vec(D_got, D_ref))
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'm=', m, ',n=', n, ',nb=', nb
        call report_case(trim(label), err, tol)
        deallocate(A0, Q, T_arr, A_ref, A_got, T_ref, T_got, D_ref, D_got)
    end do
    call report_finalize()
end program test_dorhr_col
