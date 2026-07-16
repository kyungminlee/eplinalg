! zgelsy: complex complete orthogonal factorization-based least squares.
program test_zgelsy
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_matrix_complex
    use target_lapack,   only: target_name, target_eps, target_zgelsy
    use ref_quad_lapack, only: zgelsy
    implicit none

    integer, parameter :: ms(*)  = [16, 24]
    integer, parameter :: ns(*)  = [8, 16]
    integer, parameter :: nrhs   = 2
    integer :: i, m, n, ldb, info, lwork, rank_r, rank_g
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:), work(:)
    real(ep), allocatable :: rwork(:)
    integer, allocatable :: jpvt_r(:), jpvt_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: rcond, err, tol
    character(len=48) :: label

    call report_init('zgelsy', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_complex(m, n,    A0, seed = 152101 + 47 * i)
        call gen_matrix_complex(ldb, nrhs, B0, seed = 152111 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        allocate(jpvt_r(n), jpvt_g(n), rwork(2*n))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        jpvt_r = 0; jpvt_g = 0
        rcond = -1.0_ep
        call zgelsy(m, n, nrhs, A_r, m, B_r, ldb, jpvt_r, rcond, rank_r, &
                    wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zgelsy(m, n, nrhs, A_r, m, B_r, ldb, jpvt_r, rcond, rank_r, &
                    work, lwork, rwork, info)
        deallocate(work)
        call target_zgelsy(m, n, nrhs, A_g, m, B_g, ldb, jpvt_g, rcond, rank_g, info)
        err = max_rel_err_mat_z(B_g(1:n,:), B_r(1:n,:))
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g, jpvt_r, jpvt_g, rwork)
    end do
    call report_finalize()
end program test_zgelsy
