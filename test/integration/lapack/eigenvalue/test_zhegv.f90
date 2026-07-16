program test_zhegv
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use compare,         only: max_rel_err_vec
    use test_data,       only: gen_hermitian_matrix_quad, gen_hpd_matrix_quad
    use target_lapack,   only: target_name, target_eps, target_zhegv
    use ref_quad_lapack, only: zhegv
    implicit none

    integer, parameter :: ns(*) = [8, 24, 48]
    integer :: i, n, info, lwork
    complex(ep), allocatable :: A0(:,:), B0(:,:), A_ref(:,:), B_ref(:,:), A_got(:,:), B_got(:,:), work(:)
    real(ep),    allocatable :: w_ref(:), w_got(:), rwork(:)
    complex(ep) :: wopt(1)
    real(ep) :: err, tol
    character(len=48) :: label

    call report_init('zhegv', target_name)
    do i = 1, size(ns)
        n = ns(i)
        call gen_hermitian_matrix_quad(n, A0, seed = 35001 + 47 * i)
        call gen_hpd_matrix_quad(n, B0, seed = 35011 + 47 * i)
        allocate(A_ref(n,n), B_ref(n,n), A_got(n,n), B_got(n,n))
        allocate(w_ref(n), w_got(n), rwork(max(1, 3*n - 2)))
        A_ref = A0; A_got = A0; B_ref = B0; B_got = B0
        call zhegv(1, 'N', 'U', n, A_ref, n, B_ref, n, w_ref, wopt, -1, rwork, info)
        lwork = max(1, int(real(wopt(1), ep)))
        allocate(work(lwork))
        call zhegv(1, 'N', 'U', n, A_ref, n, B_ref, n, w_ref, work, lwork, rwork, info)
        deallocate(work)
        call target_zhegv(1, 'N', 'U', n, A_got, n, B_got, n, w_got, info)
        err = max_rel_err_vec(w_got, w_ref)
        tol = 16.0_ep * real(n, ep)**2 * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(A_ref, B_ref, A_got, B_got, w_ref, w_got, rwork)
    end do
    call report_finalize()
end program test_zhegv
