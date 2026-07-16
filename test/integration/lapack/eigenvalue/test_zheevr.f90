program test_zheevr
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zheevr
    use ref_quad_lapack, only: zheevr
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    character(len=1), parameter :: jobzs(2) = ['N', 'V']
    integer :: i, n, info, jz, lwork, lrwork, liwork, mref, mgot, iwopt(1)
    complex(ep), allocatable :: A0(:,:), A_ref(:,:), A_got(:,:), z_ref(:,:), z_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    integer,     allocatable :: isuppz(:), iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), abstol, vl, vu, err, tol
    character(len=48) :: label

    call report_init('zheevr', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 33001 + 47 * i)
        abstol = 0.0_ep; vl = 0.0_ep; vu = 0.0_ep
        do jz = 1, size(jobzs)
            allocate(A_ref(n,n), A_got(n,n), w_ref(n), w_got(n))
            allocate(z_ref(n,n), z_got(n,n), isuppz(2*n))
            A_ref = A0; A_got = A0
            call zheevr(jobzs(jz), 'A', 'U', n, A_ref, n, vl, vu, 1, n, &
                        abstol, mref, w_ref, z_ref, n, isuppz, &
                        wopt, -1, rwopt, -1, iwopt, -1, info)
            lwork  = max(1, int(real(wopt(1), ep)))
            lrwork = max(1, int(rwopt(1)))
            liwork = max(1, iwopt(1))
            allocate(work(lwork), rwork(lrwork), iwork(liwork))
            call zheevr(jobzs(jz), 'A', 'U', n, A_ref, n, vl, vu, 1, n, &
                        abstol, mref, w_ref, z_ref, n, isuppz, &
                        work, lwork, rwork, lrwork, iwork, liwork, info)
            deallocate(work, rwork, iwork)
            call target_zheevr(jobzs(jz), 'A', 'U', n, A_got, n, vl, vu, 1, n, &
                               abstol, mgot, w_got, z_got, n, isuppz, info)
            err = max_rel_err_vec(w_got(1:mgot), w_ref(1:mref))
            tol = 16.0_ep * real(n, ep)**2 * target_eps
            write(label, '(a,a,a,i0)') 'jobz=', jobzs(jz), ',n=', n
            call report_case(trim(label), err, tol)
            deallocate(A_ref, A_got, w_ref, w_got, z_ref, z_got, isuppz)
        end do
    end do
    call report_finalize()
end program test_zheevr
