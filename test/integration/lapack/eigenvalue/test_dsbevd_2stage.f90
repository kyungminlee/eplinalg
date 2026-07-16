! 2-stage symmetric banded D&C eigensolver. JOBZ='N' only.
program test_dsbevd_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, pack_sym_band_quad
    use target_lapack,   only: target_name, target_eps, target_dsbevd_2stage
    use ref_quad_lapack, only: dsbevd_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer :: i, n, kd, info, ldab, lwork, liwork, iwopt(1)
    real(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:)
    real(ep), allocatable :: w_ref(:), w_got(:), z_ref(:,:), z_got(:,:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    character(len=48) :: label

    call report_init('dsbevd_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_symmetric_matrix_quad(n, A0, seed = 90001 + 47 * i)
        allocate(AB_ref(ldab, n), AB_got(ldab, n))
        allocate(w_ref(n), w_got(n), z_ref(n, n), z_got(n, n))
        call pack_sym_band_quad('U', n, kd, A0, AB_ref)
        AB_got = AB_ref
        call dsbevd_2stage('N', 'U', n, kd, AB_ref, ldab, w_ref, z_ref, n, &
                           wopt, -1, iwopt, -1, info)
        lwork  = max(1, int(wopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dsbevd_2stage('N', 'U', n, kd, AB_ref, ldab, w_ref, z_ref, n, &
                           work, lwork, iwork, liwork, info)
        deallocate(work, iwork)
        call target_dsbevd_2stage('N', 'U', n, kd, AB_got, ldab, &
                                  w_got, z_got, n, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'jobz=N,uplo=U,n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB_ref, AB_got, w_ref, w_got, z_ref, z_got)
    end do
    call report_finalize()
end program test_dsbevd_2stage
