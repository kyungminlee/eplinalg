! dgelss: SVD-based minimum-norm least squares.
program test_dgelss
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat, max_rel_err_vec
    use test_data,       only: gen_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dgelss
    use ref_quad_lapack, only: dgelss
    implicit none

    integer, parameter :: ms(*)  = [16, 24]
    integer, parameter :: ns(*)  = [8, 16]
    integer, parameter :: nrhs   = 2
    integer :: i, m, n, ldb, info, lwork, rank_r, rank_g
    real(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:)
    real(ep), allocatable :: S_r(:), S_g(:), work(:)
    real(ep) :: wopt(1), rcond, err_b, err_s, err, tol
    character(len=48) :: label

    call report_init('dgelss', target_name)
    do i = 1, size(ms)
        m = ms(i); n = ns(i); ldb = max(m, n)
        call gen_matrix_quad(m, n,    A0, seed = 151001 + 47 * i)
        call gen_matrix_quad(ldb, nrhs, B0, seed = 151011 + 47 * i)
        allocate(A_r(m,n), B_r(ldb,nrhs), A_g(m,n), B_g(ldb,nrhs))
        allocate(S_r(min(m,n)), S_g(min(m,n)))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        rcond = -1.0_ep
        call dgelss(m, n, nrhs, A_r, m, B_r, ldb, S_r, rcond, rank_r, &
                    wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgelss(m, n, nrhs, A_r, m, B_r, ldb, S_r, rcond, rank_r, &
                    work, lwork, info)
        deallocate(work)
        call target_dgelss(m, n, nrhs, A_g, m, B_g, ldb, S_g, rcond, rank_g, info)
        err_b = max_rel_err_mat(B_g(1:n,:), B_r(1:n,:))
        err_s = max_rel_err_vec(S_g, S_r)
        err = max(err_b, err_s)
        tol = 256.0_ep * real(max(m,n), ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'm=', m, ',n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, A_g, B_g, S_r, S_g)
    end do
    call report_finalize()
end program test_dgelss
