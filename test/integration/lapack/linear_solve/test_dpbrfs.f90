program test_dpbrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat
    use test_data,       only: gen_spd_matrix_quad, gen_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_dpbtrf, target_dpbtrs, target_dpbrfs
    use ref_quad_lapack, only: dpbtrf, dpbtrs, dpbrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer, parameter :: nrhs = 3
    integer :: i, n, kd, ldab, info
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A0(:,:), AB0(:,:), AFB_ref(:,:), AFB_got(:,:)
    real(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('dpbrfs', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_spd_matrix_quad(n, A0, seed = 140001 + 47 * i)
        call gen_matrix_quad(n, nrhs, B0, seed = 140011 + 47 * i)
        allocate(AB0(ldab, n), AFB_ref(ldab, n), AFB_got(ldab, n))
        call pack_sym_band_quad('U', n, kd, A0, AB0)
        AFB_ref = AB0; AFB_got = AB0
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(3*n), iwork(n), &
                 ferr(nrhs), berr(nrhs))
        call dpbtrf('U', n, kd, AFB_ref, ldab, info)
        call target_dpbtrf('U', n, kd, AFB_got, ldab, info)
        X_ref = B0; X_got = B0
        call dpbtrs('U', n, kd, nrhs, AFB_ref, ldab, X_ref, n, info)
        call target_dpbtrs('U', n, kd, nrhs, AFB_got, ldab, X_got, n, info)
        call dpbrfs('U', n, kd, nrhs, AB0, ldab, AFB_ref, ldab, B0, n, &
                    X_ref, n, ferr, berr, work, iwork, info)
        call target_dpbrfs('U', n, kd, nrhs, AB0, ldab, AFB_got, ldab, B0, n, &
                           X_got, n, ferr, berr, info)
        err = max_rel_err_mat(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB0, AFB_ref, AFB_got, X_ref, X_got, work, iwork, ferr, berr)
    end do
    call report_finalize()
end program test_dpbrfs
