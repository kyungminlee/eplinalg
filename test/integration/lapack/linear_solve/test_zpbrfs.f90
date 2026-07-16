program test_zpbrfs
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_mat_z
    use test_data,       only: gen_hpd_matrix_quad, gen_matrix_complex, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, &
                                target_zpbtrf, target_zpbtrs, target_zpbrfs
    use ref_quad_lapack, only: zpbtrf, zpbtrs, zpbrfs
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer, parameter :: nrhs = 3
    integer :: i, n, kd, ldab, info
    complex(ep), allocatable :: A0(:,:), AB0(:,:), AFB_ref(:,:), AFB_got(:,:)
    complex(ep), allocatable :: B0(:,:), X_ref(:,:), X_got(:,:), work(:)
    real(ep), allocatable :: rwork(:), ferr(:), berr(:)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zpbrfs', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hpd_matrix_quad(n, A0, seed = 141001 + 47 * i)
        call gen_matrix_complex(n, nrhs, B0, seed = 141011 + 47 * i)
        allocate(AB0(ldab, n), AFB_ref(ldab, n), AFB_got(ldab, n))
        call pack_herm_band_quad('U', n, kd, A0, AB0)
        AFB_ref = AB0; AFB_got = AB0
        allocate(X_ref(n, nrhs), X_got(n, nrhs), work(2*n), rwork(n), &
                 ferr(nrhs), berr(nrhs))
        call zpbtrf('U', n, kd, AFB_ref, ldab, info)
        call target_zpbtrf('U', n, kd, AFB_got, ldab, info)
        X_ref = B0; X_got = B0
        call zpbtrs('U', n, kd, nrhs, AFB_ref, ldab, X_ref, n, info)
        call target_zpbtrs('U', n, kd, nrhs, AFB_got, ldab, X_got, n, info)
        call zpbrfs('U', n, kd, nrhs, AB0, ldab, AFB_ref, ldab, B0, n, &
                    X_ref, n, ferr, berr, work, rwork, info)
        call target_zpbrfs('U', n, kd, nrhs, AB0, ldab, AFB_got, ldab, B0, n, &
                           X_got, n, ferr, berr, info)
        err = max_rel_err_mat_z(X_got, X_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB0, AFB_ref, AFB_got, X_ref, X_got, work, rwork, ferr, berr)
    end do
    call report_finalize()
end program test_zpbrfs
