! 2-stage symmetric banded selective eigensolver. JOBZ='N', RANGE='A'.
program test_dsbevx_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dsbevx_2stage
    use ref_quad_lapack, only: dsbevx_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer :: i, n, kd, info, ldab, lwork, mref, mgot
    real(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:)
    real(ep), allocatable :: Q_ref(:,:), Q_got(:,:)
    real(ep), allocatable :: w_ref(:), w_got(:), z_ref(:,:), z_got(:,:), work(:)
    integer,  allocatable :: iwork(:), ifail(:)
    real(ep) :: wopt(1), abstol, vl, vu, err, tol
    character(len=48) :: label

    call report_init('dsbevx_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_symmetric_matrix_quad(n, A0, seed = 92001 + 47 * i)
        abstol = 0.0_ep; vl = 0.0_ep; vu = 0.0_ep
        allocate(AB_ref(ldab, n), AB_got(ldab, n))
        allocate(Q_ref(n, n), Q_got(n, n))
        allocate(w_ref(n), w_got(n), z_ref(n, n), z_got(n, n))
        allocate(iwork(5*n), ifail(n))
        call pack_sym_band_quad('U', n, kd, A0, AB_ref)
        AB_got = AB_ref
        call dsbevx_2stage('N', 'A', 'U', n, kd, AB_ref, ldab, Q_ref, n, &
                           vl, vu, 1, n, abstol, mref, w_ref, z_ref, n, &
                           wopt, -1, iwork, ifail, info)
        lwork = max(1, int(wopt(1)))
        allocate(work(lwork))
        call dsbevx_2stage('N', 'A', 'U', n, kd, AB_ref, ldab, Q_ref, n, &
                           vl, vu, 1, n, abstol, mref, w_ref, z_ref, n, &
                           work, lwork, iwork, ifail, info)
        deallocate(work)
        call target_dsbevx_2stage('N', 'A', 'U', n, kd, AB_got, ldab, &
                                  Q_got, n, vl, vu, 1, n, abstol, mgot, &
                                  w_got, z_got, n, ifail, info)
        err = max_rel_err_vec(w_got(1:mgot), w_ref(1:mref))
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'jobz=N,range=A,uplo=U,n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB_ref, AB_got, Q_ref, Q_got, w_ref, w_got, z_ref, z_got, &
                   iwork, ifail)
    end do
    call report_finalize()
end program test_dsbevx_2stage
