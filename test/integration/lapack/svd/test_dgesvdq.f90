! dgesvdq: SVD with optional pivoted QR preconditioning.
program test_dgesvdq
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgesvdq
    use ref_quad_lapack, only: dgesvdq
    implicit none

    integer, parameter :: ms(*) = [16, 32]
    integer, parameter :: ns(*) = [8,  16]
    integer :: i, m, n, info, numrank_r, numrank_g, lwork, lrwork, liwork
    real(ep), allocatable :: A0(:,:), A_r(:,:), A_g(:,:), S_r(:), S_g(:), U(:,:), V(:,:), work(:), rwork(:)
    integer, allocatable :: iwork(:)
    ! dgesvdq workspace query writes WORK(1) and WORK(2); allocate at least 2.
    real(ep) :: wopt(2), ropt(2), err, tol
    character(len=48) :: label

    call report_init('dgesvdq', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i)
        call gen_matrix_quad(m, n, A0, seed = 25001 + 67 * i)
        allocate(A_r(m, n), A_g(m, n), S_r(min(m, n)), S_g(min(m, n)), U(m, n), V(n, n))
        liwork = m + n
        allocate(iwork(liwork))
        A_r = A0; A_g = A0
        call dgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_r, m, S_r, U, m, V, n, &
                     numrank_r, iwork, -1, wopt, -1, ropt, -1, info)
        lwork  = max(1, int(wopt(1)))
        lrwork = max(1, int(ropt(1)))
        liwork = max(liwork, iwork(1))
        deallocate(iwork); allocate(iwork(liwork), work(lwork), rwork(lrwork))
        call dgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_r, m, S_r, U, m, V, n, &
                     numrank_r, iwork, liwork, work, lwork, rwork, lrwork, info)
        deallocate(work, rwork)
        call target_dgesvdq('A', 'N', 'N', 'N', 'N', m, n, A_g, m, S_g, U, m, V, n, &
                            numrank_g, info)
        err = max_rel_err_vec(S_g, S_r)
        tol = 16.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, A_r, A_g, S_r, S_g, U, V, iwork)
    end do
    call report_finalize()
end program test_dgesvdq
