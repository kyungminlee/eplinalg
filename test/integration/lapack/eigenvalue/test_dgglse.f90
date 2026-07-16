! dgglse: equality-constrained least squares. Smoke test.
program test_dgglse
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad, gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dgglse
    use ref_quad_lapack, only: dgglse
    implicit none

    integer, parameter :: ns(*) = [12, 24, 32]
    integer :: i, n, m, p, info, lwork
    real(ep), allocatable :: A(:,:), B(:,:), c(:), d(:)
    real(ep), allocatable :: A_r(:,:), B_r(:,:), c_r(:), d_r(:), x_r(:), work(:)
    real(ep), allocatable :: A_g(:,:), B_g(:,:), c_g(:), d_g(:), x_g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dgglse', target_name)
    do i = 1, size(ns)
        n = ns(i); m = n; p = n / 2
        call gen_matrix_quad(m, n, A, seed = 143001 + 47 * i)
        call gen_matrix_quad(p, n, B, seed = 143011 + 47 * i)
        allocate(c(m), d(p))
        call gen_vector_quad(m, c, seed = 143021 + 47 * i)
        call gen_vector_quad(p, d, seed = 143031 + 47 * i)
        allocate(A_r(m,n), B_r(p,n), c_r(m), d_r(p), x_r(n))
        allocate(A_g(m,n), B_g(p,n), c_g(m), d_g(p), x_g(n))
        A_r = A; B_r = B; c_r = c; d_r = d
        A_g = A; B_g = B; c_g = c; d_g = d
        call dgglse(m, n, p, A_r, m, B_r, p, c_r, d_r, x_r, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgglse(m, n, p, A_r, m, B_r, p, c_r, d_r, x_r, work, lwork, info)
        deallocate(work)
        call target_dgglse(m, n, p, A_g, m, B_g, p, c_g, d_g, x_g, info)
        err = max_rel_err_vec(x_g, x_r)
        tol = 64.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, c, d, A_r, B_r, c_r, d_r, x_r, A_g, B_g, c_g, d_g, x_g)
    end do
    call report_finalize()
end program test_dgglse
