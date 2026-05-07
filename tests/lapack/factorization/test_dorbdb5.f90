! dorbdb5: project [X1; X2] onto the column orthogonal complement of [Q1; Q2].
program test_dorbdb5
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_matrix_quad, gen_vector_quad
    use target_lapack,   only: target_name, target_eps, target_dorbdb5
    use ref_quad_lapack, only: dgeqrf, dorgqr, dorbdb5
    implicit none

    integer, parameter :: ms(*) = [12, 16]
    integer :: i, m, m1, m2, n, info, lwork
    real(ep), allocatable :: A(:,:), tau(:), work(:)
    real(ep), allocatable :: Q1(:,:), Q2(:,:), X1r(:), X2r(:), X1g(:), X2g(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dorbdb5', target_name)
    do i = 1, size(ms)
        m = ms(i); m1 = m/2; m2 = m - m1; n = m/3
        ! Build orthogonal m×m, slice first n columns into [Q1; Q2].
        call gen_matrix_quad(m, m, A, seed = 24701 + 79 * i)
        allocate(tau(m))
        call dgeqrf(m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dgeqrf(m, m, A, m, tau, work, lwork, info)
        deallocate(work)
        call dorgqr(m, m, m, A, m, tau, wopt, -1, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dorgqr(m, m, m, A, m, tau, work, lwork, info)
        deallocate(work, tau)
        allocate(Q1(m1, n), Q2(m2, n))
        Q1 = A(1:m1, 1:n); Q2 = A(m1+1:m, 1:n)
        deallocate(A)
        allocate(X1r(m1), X2r(m2), X1g(m1), X2g(m2))
        call gen_vector_quad(m1, X1r, seed = 24711 + 79 * i)
        call gen_vector_quad(m2, X2r, seed = 24721 + 79 * i)
        X1g = X1r; X2g = X2r
        ! dorbdb5 has no workspace query; LWORK >= N is the spec.
        lwork = n
        allocate(work(lwork))
        call dorbdb5(m1, m2, n, X1r, 1, X2r, 1, Q1, m1, Q2, m2, work, lwork, info)
        deallocate(work)
        call target_dorbdb5(m1, m2, n, X1g, 1, X2g, 1, Q1, m1, Q2, m2, info)
        err = max(max_rel_err_vec(X1g, X1r), max_rel_err_vec(X2g, X2r))
        tol = 16.0_ep * real(m, ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0,a,i0)') 'm1=', m1, ',m2=', m2, ',n=', n, ',m=', m
        call report_case(trim(label), err, tol)
        deallocate(Q1, Q2, X1r, X2r, X1g, X2g)
    end do
    call report_finalize()
end program test_dorbdb5
