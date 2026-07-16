! 2-stage Hermitian banded eigensolver. JOBZ='N' only. Adds an LWORK
! parameter not present in the non-2-stage variant.
program test_zhbev_2stage
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zhbev_2stage
    use ref_quad_lapack, only: zhbev_2stage
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    integer :: i, n, kd, info, ldab, lwork
    complex(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:)
    complex(ep), allocatable :: z_ref(:,:), z_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhbev_2stage', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hermitian_matrix_quad(n, A0, seed = 89001 + 47 * i)
        allocate(AB_ref(ldab, n), AB_got(ldab, n))
        allocate(w_ref(n), w_got(n), z_ref(n, n), z_got(n, n))
        allocate(rwork(max(1, n)))
        call pack_herm_band_quad('U', n, kd, A0, AB_ref)
        AB_got = AB_ref
        call zhbev_2stage('N', 'U', n, kd, AB_ref, ldab, &
                          w_ref, z_ref, n, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhbev_2stage('N', 'U', n, kd, AB_ref, ldab, &
                          w_ref, z_ref, n, work, lwork, rwork, info)
        deallocate(work)
        call target_zhbev_2stage('N', 'U', n, kd, AB_got, ldab, &
                                 w_got, z_got, n, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0,a,i0)') 'jobz=N,uplo=U,n=', n, ',kd=', kd
        call report_case(trim(label), err, tol)
        deallocate(AB_ref, AB_got, w_ref, w_got, z_ref, z_got, rwork)
    end do
    call report_finalize()
end program test_zhbev_2stage
