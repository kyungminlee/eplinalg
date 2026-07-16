program test_zhpgvd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, gen_hpd_matrix_quad, &
                                 pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhpgvd
    use ref_quad_lapack, only: zhpgvd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, np, info, lwork, lrwork, liwork
    complex(ep), allocatable :: A(:,:), B(:,:), work(:)
    complex(ep), allocatable :: AP_ref(:), BP_ref(:), AP_got(:), BP_got(:), Z(:,:)
    real(ep), allocatable :: W_ref(:), W_got(:), rwork(:)
    integer,  allocatable :: iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    integer  :: iwopt(1)
    character(len=48) :: label

    call report_init('zhpgvd', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A, seed = 260061 + 47 * i)
        call gen_hpd_matrix_quad(n, B, seed = 260071 + 47 * i)
        allocate(AP_ref(np), BP_ref(np), AP_got(np), BP_got(np))
        call pack_herm_packed_quad('U', n, A, AP_ref); AP_got = AP_ref
        call pack_herm_packed_quad('U', n, B, BP_ref); BP_got = BP_ref
        allocate(W_ref(n), W_got(n), Z(1, 1))
        call zhpgvd(1, 'N', 'U', n, AP_ref, BP_ref, W_ref, Z, 1, wopt, -1, &
                    rwopt, -1, iwopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep))); lrwork = max(1, int(rwopt(1)))
        liwork = max(1, iwopt(1))
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        call zhpgvd(1, 'N', 'U', n, AP_ref, BP_ref, W_ref, Z, 1, work, lwork, &
                    rwork, lrwork, iwork, liwork, info)
        call target_zhpgvd(1, 'N', 'U', n, AP_got, BP_got, W_got, Z, 1, info)
        err = max_rel_err_vec(W_got, W_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, B, AP_ref, BP_ref, AP_got, BP_got, W_ref, W_got, Z, work, rwork, iwork)
    end do
    call report_finalize()
end program test_zhpgvd
