program test_zhbevd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_band_quad
    use target_lapack,   only: target_name, target_eps, target_zhbevd
    use ref_quad_lapack, only: zhbevd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 64]
    integer, parameter :: kds(*) = [2, 4, 8]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    integer :: i, n, kd, info, jz, ldab, lwork, lrwork, liwork, iwopt(1)
    integer, allocatable :: iwork(:)
    complex(ep), allocatable :: A0(:,:), AB_ref(:,:), AB_got(:,:)
    complex(ep), allocatable :: z_ref(:,:), z_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    character(len=48) :: label

    call report_init('zhbevd', target_name)
    do i = 1, size(ns)
        n = ns(i); kd = kds(i); ldab = kd + 1
        call gen_hermitian_matrix_quad(n, A0, seed = 51001 + 47 * i)
        do jz = 1, size(jobzs)
            allocate(AB_ref(ldab, n), AB_got(ldab, n))
            allocate(w_ref(n), w_got(n), z_ref(n, n), z_got(n, n))
            call pack_herm_band_quad('U', n, kd, A0, AB_ref)
            AB_got = AB_ref
            call zhbevd(jobzs(jz), 'U', n, kd, AB_ref, ldab, w_ref, z_ref, n, &
                        wopt, -1, rwopt, -1, iwopt, -1, info)
            lwork  = max(1, int(real(wopt(1), ep)))
            lrwork = max(1, int(rwopt(1)))
            liwork = max(1, iwopt(1))
            allocate(work(lwork), rwork(lrwork), iwork(liwork))
            call zhbevd(jobzs(jz), 'U', n, kd, AB_ref, ldab, w_ref, z_ref, n, &
                        work, lwork, rwork, lrwork, iwork, liwork, info)
            deallocate(work, rwork, iwork)
            call target_zhbevd(jobzs(jz), 'U', n, kd, AB_got, ldab, &
                               w_got, z_got, n, info)
            err = max_rel_err_vec(w_got, w_ref)
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0,a,i0)') 'jobz=', jobzs(jz), ',n=', n, ',kd=', kd
            call report_case(trim(label), err, tol)
            deallocate(AB_ref, AB_got, w_ref, w_got, z_ref, z_got)
        end do
    end do
    call report_finalize()
end program test_zhbevd
