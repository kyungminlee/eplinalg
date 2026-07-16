! zhegvx: complex Hermitian gen-eig expert driver.
program test_zhegvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhegvx
    use ref_quad_lapack, only: zhegvx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info, m_r, m_g, lwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_r(:,:), B_r(:,:), A_g(:,:), B_g(:,:)
    complex(ep), allocatable :: Z_r(:,:), Z_g(:,:), work(:)
    real(ep), allocatable :: w_r(:), w_g(:), rwork(:)
    integer, allocatable :: iwork(:), ifail_r(:), ifail_g(:)
    complex(ep) :: wopt(1)
    real(ep) :: vl, vu, abstol, err, tol
    character(len=48) :: label

    call report_init('zhegvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 22501 + 47 * i)
        call gen_hpd_matrix_quad(n, B0, seed = 22511 + 47 * i)
        allocate(A_r(n,n), B_r(n,n), w_r(n), Z_r(n,n), ifail_r(n))
        allocate(A_g(n,n), B_g(n,n), w_g(n), Z_g(n,n), ifail_g(n))
        A_r = A0; B_r = B0; A_g = A0; B_g = B0
        vl = 0.0_ep; vu = 0.0_ep; abstol = 0.0_ep
        allocate(rwork(7*n), iwork(5*n))
        call zhegvx(1, 'N', 'A', 'U', n, A_r, n, B_r, n, vl, vu, 1, n, abstol, &
                    m_r, w_r, Z_r, n, wopt, -1, rwork, iwork, ifail_r, info)
        lwork = max(1, int(real(wopt(1), ep))); allocate(work(lwork))
        call zhegvx(1, 'N', 'A', 'U', n, A_r, n, B_r, n, vl, vu, 1, n, abstol, &
                    m_r, w_r, Z_r, n, work, lwork, rwork, iwork, ifail_r, info)
        deallocate(work, rwork, iwork)
        call target_zhegvx(1, 'N', 'A', 'U', n, A_g, n, B_g, n, vl, vu, 1, n, abstol, &
                           m_g, w_g, Z_g, n, ifail_g, info)
        err = max_rel_err_vec(w_g(1:m_g), w_r(1:m_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_r, B_r, w_r, Z_r, ifail_r, A_g, B_g, w_g, Z_g, ifail_g)
    end do
    call report_finalize()
end program test_zhegvx
