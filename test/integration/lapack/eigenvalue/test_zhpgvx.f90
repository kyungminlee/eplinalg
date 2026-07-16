! zhpgvx: packed complex Hermitian gen-eig expert driver.
program test_zhpgvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hpd_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhpgvx
    use ref_quad_lapack, only: zhpgvx
    implicit none

    integer, parameter :: ns(*) = [12, 24, 48]
    integer :: i, n, info, m_r, m_g
    complex(ep), allocatable :: A0(:,:), B0(:,:), AP0(:), BP0(:)
    complex(ep), allocatable :: AP_r(:), BP_r(:), Z_r(:,:), work(:)
    complex(ep), allocatable :: AP_g(:), BP_g(:), Z_g(:,:)
    real(ep), allocatable :: w_r(:), w_g(:), rwork(:)
    integer, allocatable :: iwork(:), ifail_r(:), ifail_g(:)
    real(ep) :: vl, vu, abstol, err, tol
    character(len=48) :: label

    call report_init('zhpgvx', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hpd_matrix_quad(n, A0, seed = 22301 + 47 * i)
        call gen_hpd_matrix_quad(n, B0, seed = 22311 + 47 * i)
        allocate(AP0(n*(n+1)/2), BP0(n*(n+1)/2))
        call pack_herm_packed_quad('U', n, A0, AP0)
        call pack_herm_packed_quad('U', n, B0, BP0)
        allocate(AP_r(n*(n+1)/2), BP_r(n*(n+1)/2), w_r(n), Z_r(n,n), ifail_r(n))
        allocate(AP_g(n*(n+1)/2), BP_g(n*(n+1)/2), w_g(n), Z_g(n,n), ifail_g(n))
        AP_r = AP0; BP_r = BP0; AP_g = AP0; BP_g = BP0
        vl = 0.0_ep; vu = 0.0_ep; abstol = 0.0_ep
        allocate(work(2*n), rwork(7*n), iwork(5*n))
        call zhpgvx(1, 'N', 'A', 'U', n, AP_r, BP_r, vl, vu, 1, n, abstol, m_r, &
                    w_r, Z_r, n, work, rwork, iwork, ifail_r, info)
        deallocate(work, rwork, iwork)
        call target_zhpgvx(1, 'N', 'A', 'U', n, AP_g, BP_g, vl, vu, 1, n, abstol, &
                           m_g, w_g, Z_g, n, ifail_g, info)
        err = max_rel_err_vec(w_g(1:m_g), w_r(1:m_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, AP0, BP0, AP_r, BP_r, w_r, Z_r, ifail_r, &
                   AP_g, BP_g, w_g, Z_g, ifail_g)
    end do
    call report_finalize()
end program test_zhpgvx
