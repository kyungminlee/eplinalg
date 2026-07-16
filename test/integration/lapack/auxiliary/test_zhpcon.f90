program test_zhpcon
    use prec_kinds,      only: ep
    use prec_report,     only: report_init, report_case, report_finalize
    use test_data,       only: gen_hermitian_matrix_quad, pack_herm_packed_quad
    use target_lapack,   only: target_name, target_eps, target_zhptrf, target_zhpcon
    use ref_quad_lapack, only: zhptrf, zhpcon
    implicit none
    integer, parameter :: ns(*) = [16, 32, 64]
    integer :: i, n, info, np, j
    integer, allocatable :: ipiv_ref(:), ipiv_got(:)
    complex(ep), allocatable :: A0(:,:), AP_ref(:), AP_got(:), work(:)
    real(ep) :: anorm, rcond_ref, rcond_got, err, tol
    character(len=48) :: label
    call report_init('zhpcon', target_name)
    do i = 1, size(ns)
        n = ns(i); np = n*(n+1)/2
        call gen_hermitian_matrix_quad(n, A0, seed = 73001 + 47 * i)
        do j = 1, n; A0(j,j) = A0(j,j) + cmplx(real(n, ep), 0.0_ep, ep); end do
        anorm = maxval(sum(abs(A0), dim=1))
        allocate(AP_ref(np), AP_got(np), ipiv_ref(n), ipiv_got(n), work(2*n))
        call pack_herm_packed_quad('U', n, A0, AP_ref); AP_got = AP_ref
        call zhptrf('U', n, AP_ref, ipiv_ref, info)
        call target_zhptrf('U', n, AP_got, ipiv_got, info)
        call zhpcon('U', n, AP_ref, ipiv_ref, anorm, rcond_ref, work, info)
        call target_zhpcon('U', n, AP_got, ipiv_got, anorm, rcond_got, info)
        err = abs(rcond_got - rcond_ref) / max(abs(rcond_ref), tiny(1.0_ep))
        tol = 100.0_ep * real(n, ep) * target_eps
        write(label, '(a,i0)') 'n=', n
        call report_case(trim(label), err, tol)
        deallocate(AP_ref, AP_got, ipiv_ref, ipiv_got, work)
    end do
    call report_finalize()
end program test_zhpcon
