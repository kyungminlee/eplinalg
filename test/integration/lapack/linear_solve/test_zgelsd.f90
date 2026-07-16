! zgelsd: complex divide-and-conquer SVD-based least squares.
program test_zgelsd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z, max_rel_err_vec
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgelsd
    use ref_quad_lapack, only: zgelsd
    implicit none

    integer, parameter :: ms(*)  = [16, 24]
    integer, parameter :: ns(*)  = [8, 16]
    integer, parameter :: nrhs   = 2
    integer :: i, m, n, ldb, info, lwork, lrwork, rank_r, rank_g
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:), work(:)
    real(ep), allocatable :: S_r(:), S_g(:), rwork(:)
    integer, allocatable :: iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), rcond, err_b, err_s, err, tol
    integer :: iwopt(1)
    character(len=48) :: label

    call report_init('zgelsd', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_complex(m, n,    A0, seed = 150101 + 47 * i)
        call gen_matrix_complex(ldb, nrhs, B0, seed = 150111 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        allocate(S_r(min(m,n)), S_g(min(m,n)))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        rcond = -1.0_ep
        call zgelsd(m, n, nrhs, A_r, m, B_r, ldb, S_r, rcond, rank_r, &
                    wopt, -1, rwopt, iwopt, info)
        lwork = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1)))
        allocate(work(lwork), rwork(lrwork), iwork(max(1, iwopt(1))))
        call zgelsd(m, n, nrhs, A_r, m, B_r, ldb, S_r, rcond, rank_r, &
                    work, lwork, rwork, iwork, info)
        deallocate(work, rwork, iwork)
        call target_zgelsd(m, n, nrhs, A_g, m, B_g, ldb, S_g, rcond, rank_g, info)
        err_b = max_rel_err_mat_z(B_g(1:n,:), B_r(1:n,:))
        err_s = max_rel_err_vec(S_g, S_r)
        err = max(err_b, err_s)
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g, S_r, S_g)
    end do
    call report_finalize()
end program test_zgelsd
