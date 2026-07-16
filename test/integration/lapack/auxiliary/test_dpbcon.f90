program test_dpbcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_spd_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dpbtrf, target_dpbcon
    use ref_quad_lapack, only: dpbtrf, dpbcon
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer :: i, n, kd, ldab, info
    integer, allocatable :: iwork(:)
    real(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label

    call report_init('dpbcon', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_spd_matrix_quad(n, A0, seed = 120001 + 47 * i)
        allocate(AB_ref(ldab, n), AB_got(ldab, n), work(3*n), iwork(n))
        call pack_sym_band_quad('U', n, kd, A0, AB_ref); AB_got = AB_ref
        anorm = maxval(sum(abs(AB_ref), dim=1))
        call dpbtrf('U', n, kd, AB_ref, ldab, info)
        call target_dpbtrf('U', n, kd, AB_got, ldab, info)
        call dpbcon('U', n, kd, AB_ref, ldab, anorm, rcond_ref, work, iwork, info)
        call target_dpbcon('U', n, kd, AB_got, ldab, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 1000.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0,a,i0)') 'n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB_ref, AB_got, work, iwork)
    end do
    call report_finalize()
end program test_dpbcon
