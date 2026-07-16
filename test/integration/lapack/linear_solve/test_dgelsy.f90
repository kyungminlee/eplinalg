! dgelsy: complete orthogonal factorization-based least squares.
program test_dgelsy
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgelsy
    use ref_quad_lapack, only: dgelsy
    implicit none

    integer, parameter :: ms(*)  = [16, 24]
    integer, parameter :: ns(*)  = [8, 16]
    integer, parameter :: nrhs   = 2
    integer :: i, m, n, ldb, info, lwork, rank_r, rank_g
    real(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:), work(:)
    integer, allocatable :: jpvt_r(:), jpvt_g(:)
    real(ep) :: wopt(1), rcond, err, tol
    character(len=48) :: label

    call report_init('dgelsy', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_quad(m, n,    A0, seed = 152001 + 47 * i)
        call gen_matrix_quad(ldb, nrhs, B0, seed = 152011 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        allocate(jpvt_r(n), jpvt_g(n))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        jpvt_r = 0; jpvt_g = 0
        rcond = -1.0_ep
        call dgelsy(m, n, nrhs, A_r, m, B_r, ldb, jpvt_r, rcond, rank_r, &
                    wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgelsy(m, n, nrhs, A_r, m, B_r, ldb, jpvt_r, rcond, rank_r, &
                    work, lwork, info)
        deallocate(work)
        call target_dgelsy(m, n, nrhs, A_g, m, B_g, ldb, jpvt_g, rcond, rank_g, info)
        err = max_rel_err_mat(B_g(1:n,:), B_r(1:n,:))
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g, jpvt_r, jpvt_g)
    end do
    call report_finalize()
end program test_dgelsy
