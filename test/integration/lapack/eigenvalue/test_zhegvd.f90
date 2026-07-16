program test_zhegvd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhegvd
    use ref_quad_lapack, only: zhegvd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork, lrwork, liwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    complex(ep), allocatable :: work(:)
    real(ep), allocatable :: W_ref(:), W_got(:), rwork(:)
    integer,  allocatable :: iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    integer  :: iwopt(1)
    character(len=48) :: label

    call report_init('zhegvd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 260021 + 47 * i)
        call gen_hpd_matrix_quad(n, B0, seed = 260031 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n), W_ref(n), W_got(n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call zhegvd(1, 'N', 'U', n, A_ref, n, B_ref, n, W_ref, wopt, -1, rwopt, -1, &
                    iwopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); lrwork = max(1, int(rwopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        call zhegvd(1, 'N', 'U', n, A_ref, n, B_ref, n, W_ref, work, lwork, rwork, lrwork, &
                    iwork, liwork, info)
        call target_zhegvd(1, 'N', 'U', n, A_got, n, B_got, n, W_got, info)
        err = max_rel_err_vec(W_got, W_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, W_ref, W_got, work, rwork, iwork)
    end do
    call report_finalize()
end program test_zhegvd
