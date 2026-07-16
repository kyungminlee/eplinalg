! dsbgvx: real symmetric banded gen-eig expert driver.
program test_dsbgvx
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_spd_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dsbgvx
    use ref_quad_lapack, only: dsbgvx
    implicit none

    integer, parameter :: ns(*)  = [12, 24]
    integer, parameter :: kas(*) = [3, 5]
    integer, parameter :: kbs(*) = [3, 5]
    integer :: i, n, ka, kb, ldab, ldbb, info, m_r, m_g
    real(ep), allocatable :: A0(:,:), B0(:,:), AB0(:,:), BB0(:,:)
    real(ep), allocatable :: AB_r(:,:), BB_r(:,:), Q_r(:,:), w_r(:), Z_r(:,:)
    real(ep), allocatable :: AB_g(:,:), BB_g(:,:), Q_g(:,:), w_g(:), Z_g(:,:)
    real(ep), allocatable :: work(:)
    integer, allocatable :: iwork(:), ifail_r(:), ifail_g(:)
    real(ep) :: vl, vu, abstol, err, tol
    character(len=48) :: label

    call report_init('dsbgvx', target_name)
    do i = 1, size(ns)
        n = ns(i); ka = kas(i); kb = kbs(i)
        ldab = ka + 1; ldbb = kb + 1
        call gen_spd_matrix_quad(n, A0, seed = 22001 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 22011 + 47 * i)
        allocate(AB0(ldab,n), BB0(ldbb,n))
        call pack_sym_band_quad('U', n, ka, A0, AB0)
        call pack_sym_band_quad('U', n, kb, B0, BB0)
        allocate(AB_r(ldab,n), BB_r(ldbb,n), Q_r(n,n), w_r(n), Z_r(n,n), ifail_r(n))
        allocate(AB_g(ldab,n), BB_g(ldbb,n), Q_g(n,n), w_g(n), Z_g(n,n), ifail_g(n))
        AB_r = AB0; BB_r = BB0; AB_g = AB0; BB_g = BB0
        vl = 0.0_ep; vu = 0.0_ep; abstol = 0.0_ep
        allocate(work(7*n), iwork(5*n))
        call dsbgvx('N', 'A', 'U', n, ka, kb, AB_r, ldab, BB_r, ldbb, &
                    Q_r, n, vl, vu, 1, n, abstol, m_r, w_r, Z_r, n, &
                    work, iwork, ifail_r, info)
        deallocate(work, iwork)
        call target_dsbgvx('N', 'A', 'U', n, ka, kb, AB_g, ldab, BB_g, ldbb, &
                           Q_g, n, vl, vu, 1, n, abstol, m_g, w_g, Z_g, n, &
                           ifail_g, info)
        err = max_rel_err_vec(w_g(1:m_g), w_r(1:m_r))
        tol = 256.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0,a,i0)') 'n=', n, ',ka=', ka, ',kb=', kb
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, AB0, BB0, AB_r, BB_r, Q_r, w_r, Z_r, ifail_r, &
                   AB_g, BB_g, Q_g, w_g, Z_g, ifail_g)
    end do
    call report_finalize()
end program test_dsbgvx
