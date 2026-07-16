program test_zstedc
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zstedc
    use ref_quad_lapack, only: zhetrd, zstedc
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork, lrwork, liwork
    complex(ep), allocatable :: A(:,:), tau(:), work(:), Z(:,:)
    real(ep), allocatable :: D0(:), E0(:), rwork(:)
    real(ep), allocatable :: D_ref(:), E_ref(:), D_got(:), E_got(:)
    integer,  allocatable :: iwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: rwopt(1), err, tol
    integer  :: iwopt(1)
    character(len=48) :: label

    call report_init('zstedc', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A, seed = 280041 + 47 * i)
        allocate(D0(n), E0(n-1), tau(n-1))
        call zhetrd('U', n, A, n, D0, E0, tau, wopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhetrd('U', n, A, n, D0, E0, tau, work, lwork, info)
        deallocate(work)
        allocate(D_ref(n), E_ref(n-1), D_got(n), E_got(n-1), Z(1, 1))
        D_ref = D0; E_ref = E0; D_got = D0; E_got = E0
        call zstedc('N', n, D_ref, E_ref, Z, 1, wopt, -1, rwopt, -1, iwopt, -1, info)
        lwork = max(1, int(real(wopt(1), ep)))
        lrwork = max(1, int(rwopt(1))); liwork = max(1, iwopt(1))
        allocate(work(lwork), rwork(lrwork), iwork(liwork))
        call zstedc('N', n, D_ref, E_ref, Z, 1, work, lwork, rwork, lrwork, &
                    iwork, liwork, info)
        call target_zstedc('N', n, D_got, E_got, Z, 1, info)
        err = max_rel_err_vec(D_got, D_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A, D0, E0, tau, D_ref, E_ref, D_got, E_got, Z, work, rwork, iwork)
    end do
    call report_finalize()
end program test_zstedc
