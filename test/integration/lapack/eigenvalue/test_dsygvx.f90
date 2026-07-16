! dsygvx: real symmetric gen-eig expert driver.
program test_dsygvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsygvx
    use ref_quad_lapack, only: dsygvx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info, m_r, m_g, lwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:)
    real(ep), allocatable :: w_r(:), Z_r(:,:), w_g(:), Z_g(:,:), work(:)
    integer, allocatable :: iwork(:), ifail_r(:), ifail_g(:)
    real(ep) :: vl, vu, abstol, wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsygvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_spd_matrix_quad(n, A0, seed = 22401 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 22411 + 47 * i)
        allocate(A_r(n,n), B_r(n,n), w_r(n), Z_r(n,n), ifail_r(n))
        allocate(A_g(n,n), B_g(n,n), w_g(n), Z_g(n,n), ifail_g(n))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        vl = 0.0_ep; vu = 0.0_ep; abstol = 0.0_ep
        allocate(iwork(5*n))
        call dsygvx(1, 'N', 'A', 'U', n, A_r, n, B_r, n, vl, vu, 1, n, abstol, &
                    m_r, w_r, Z_r, n, wopt, -1, iwork, ifail_r, info)
        lwork = max(1, int(wopt(1))); allocate(work(lwork))
        call dsygvx(1, 'N', 'A', 'U', n, A_r, n, B_r, n, vl, vu, 1, n, abstol, &
                    m_r, w_r, Z_r, n, work, lwork, iwork, ifail_r, info)
        deallocate(work, iwork)
        call target_dsygvx(1, 'N', 'A', 'U', n, A_g, n, B_g, n, vl, vu, 1, n, abstol, &
                           m_g, w_g, Z_g, n, ifail_g, info)
        err = max_rel_err_vec(w_g(1:m_g), w_r(1:m_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, w_r, Z_r, ifail_r, A_g, B_g, w_g, Z_g, ifail_g)
    end do
    call report_finalize()
end program test_dsygvx
