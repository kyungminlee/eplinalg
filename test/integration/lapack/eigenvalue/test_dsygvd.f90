program test_dsygvd
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_symmetric_matrix_quad, gen_spd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_dsygvd
    use ref_quad_lapack, only: dsygvd
    implicit none

    integer, parameter :: ns(*) = [16, 32, 48]
    integer :: i, n, info, lwork, liwork
    real(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:)
    real(ep), allocatable :: W_ref(:), W_got(:), work(:)
    integer,  allocatable :: iwork(:)
    real(ep) :: wopt(1), err, tol
    integer  :: iwopt(1)
    character(len=48) :: label

    call report_init('dsygvd', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_symmetric_matrix_quad(n, A0, seed = 260001 + 47 * i)
        call gen_spd_matrix_quad(n, B0, seed = 260011 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n), W_ref(n), W_got(n))
        A_ref = A0; B_ref = B0; A_got = A0; B_got = B0
        call dsygvd(1, 'N', 'U', n, A_ref, n, B_ref, n, W_ref, wopt, -1, iwopt, -1, info)
        lwork = max(1, int(wopt(1))); liwork = max(1, iwopt(1))
        allocate(work(lwork), iwork(liwork))
        call dsygvd(1, 'N', 'U', n, A_ref, n, B_ref, n, W_ref, work, lwork, iwork, liwork, info)
        call target_dsygvd(1, 'N', 'U', n, A_got, n, B_got, n, W_got, info)
        err = max_rel_err_vec(W_got, W_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A0, B0, A_ref, B_ref, A_got, B_got, W_ref, W_got, work, iwork)
    end do
    call report_finalize()
end program test_dsygvd
